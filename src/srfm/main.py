"""This code defines one function to run srfm.

    Note that the SRFM can be run interactively (i.e. used as a package, use the
    inside of run_srfm() as and example), or with a driver table.

- Name: main
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 25 Nov 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from . import utilities
from . import forward_model
from . import rfm_functions
from . import layer
import os
import datetime
import warnings
from . import rfm_helper
from .input_schema import validate_srfm_inputs
from netCDF4 import Dataset
import json
import copy
from collections.abc import Mapping
from numbers import Real
from .spectral_fields import SpectralField


def _resolve_output_geometry(output_format, requested, include_toa, rfm_output):
    """Resolve requested output values and altitude rows before DISORT setup.

    Args:
        output_format (str): ``"altitude"`` or direct optical-depth ``"tau"``.
        requested (array-like): Requested altitude or optical-depth values.
        include_toa (bool): Append the top-of-atmosphere output when true.
        rfm_output (OpticalDepthGrid | pandas.DataFrame): Atmospheric profile.

    Returns:
        tuple[numpy.ndarray, list[int | None] | None]: Resolved values and the
        corresponding altitude-row indices.

    Raises:
        ValueError: If an altitude lies outside the atmospheric grid.
    """
    output_values = list(requested)
    altitude_rows = None

    if output_format == "altitude":
        if hasattr(rfm_output, "altitude_lower"):
            lower_altitudes = np.asarray(rfm_output.altitude_lower, dtype=float)
            upper_altitudes = np.asarray(rfm_output.altitude_upper, dtype=float)
        else:
            lower_altitudes = rfm_output["h_lower (km)"].to_numpy(dtype=float)
            upper_altitudes = rfm_output["h_upper (km)"].to_numpy(dtype=float)
        toa_altitude = float(upper_altitudes.max())
        bottom_altitude = float(lower_altitudes.min())
        if any(
            altitude < bottom_altitude or altitude > toa_altitude
            for altitude in output_values
        ):
            raise ValueError(
                "Requested output altitudes must lie within the atmospheric grid "
                f"({bottom_altitude:g} to {toa_altitude:g} km)."
            )
        altitude_rows = []
        matched_altitudes = []
        for altitude in output_values:
            if np.isclose(altitude, toa_altitude):
                altitude_rows.append(None)
                matched_altitudes.append(toa_altitude)
            else:
                row = int(np.abs(lower_altitudes - altitude).argmin())
                altitude_rows.append(row)
                matched_altitudes.append(float(lower_altitudes[row]))
        output_values = matched_altitudes
        if include_toa and None not in altitude_rows:
            altitude_rows.append(None)
            output_values.append(toa_altitude)
    elif include_toa and not np.any(np.isclose(output_values, 0.0)):
        output_values.append(0.0)

    return np.asarray(output_values, dtype=float), altitude_rows


def _resolve_retained_outputs(values):
    """Derive the requested and temporary SRFM output arrays.

    Args:
        values (Mapping): Validated SRFM input values.

    Returns:
        tuple[set[str], set[str]]: Publicly retained semantic names and raw DISORT
        arrays required while the calculation is running.

    Raises:
        ValueError: If the retention collection contains an unsupported name.
    """
    aliases = {"rad": "uu", "radiance": "uu"}
    raw_names = {
        "rfldir",
        "rfldn",
        "flup",
        "dfdt",
        "uavg",
        "uu",
        "albmed",
        "trnmed",
    }
    configured = values["retain_outputs"]
    requested = {aliases.get(name, name) for name in configured}
    unknown = requested - raw_names - {"bbt"}
    if unknown:
        raise ValueError("Unknown retained output name(s): " + ", ".join(sorted(unknown)))
    runtime = requested & raw_names
    if "bbt" in requested or values.get("convolve_iasi"):
        runtime.add("uu")
    return requested, runtime


def _interpolate_particle_block(particle_layers, wavenumber_cm_inverse, nmom):
    """Interpolate every particle layer for one bounded spectral block.

    Args:
        particle_layers (Mapping[str, Layer]): Mie and prescribed particle layers.
        wavenumber_cm_inverse (array-like): Wavenumbers for the current block.
        nmom (int): Highest normalized Legendre-moment order requested by DISORT.

    Returns:
        dict[str, dict[str, numpy.ndarray]]: Column optical depth,
        single-scattering albedo, and normalized Legendre coefficients for each
        named layer.
    """
    return {
        layer_name: particle_layer.column_optical_properties(
            wavenumber_cm_inverse, nmom=nmom
        )
        for layer_name, particle_layer in particle_layers.items()
    }


def _prepare_boundary_spectral_fields(values, computational_wavenumber_cm_inverse):
    """Prepare scalar/spectral albedo and an optional custom solar spectrum.

    Both fields are loaded and coverage-checked before any RFM or DISORT side
    effect. Custom solar values are converted once to beam-normal
    ``W m-2 (cm-1)-1`` at their source grid.

    Args:
        values (Mapping): Validated runner inputs.
        computational_wavenumber_cm_inverse (array-like): Complete requested
            computational grid in ``cm-1``.

    Returns:
        tuple[float | SpectralField, SpectralField | None]: Prepared albedo and
        optional custom solar field.
    """
    albedo_input = values["albedo"]
    if isinstance(albedo_input, Real) and not isinstance(albedo_input, (bool, np.bool_)):
        prepared_albedo = float(albedo_input)
    else:
        prepared_albedo = SpectralField.from_specification(
            albedo_input,
            "albedo",
            minimum=0,
            maximum=1,
        )
        prepared_albedo.validate_coverage(computational_wavenumber_cm_inverse)

    solar_input = values.get("solar_spectrum")
    prepared_solar = None
    if solar_input is not None:
        prepared_solar = SpectralField.from_specification(
            solar_input,
            "solar_spectrum",
            minimum=0,
            value_units=solar_input.get("value_units"),
            solar_density=True,
        )
        prepared_solar.validate_coverage(computational_wavenumber_cm_inverse)
    return prepared_albedo, prepared_solar


def _prepare_solar_spectral_irradiance(
    prepared_custom_solar,
    computational_wavenumber_cm_inverse,
    year_day,
    *,
    sun,
):
    """Select the custom or unchanged built-in DISORT FBEAM spectrum.

    Args:
        prepared_custom_solar (SpectralField | None): Optional validated custom
            beam-normal spectrum.
        computational_wavenumber_cm_inverse (array-like): Actual RFM grid.
        year_day (int): Day of year used only by the built-in Gueymard pathway.
        sun (bool): Whether the direct solar beam is enabled.

    Returns:
        numpy.ndarray: Beam-normal spectral irradiance in
        ``W m-2 (cm-1)-1``. Disabled sunlight returns zeros.
    """
    wavenumber = np.asarray(computational_wavenumber_cm_inverse, dtype=float)
    if not sun:
        return np.zeros(wavenumber.shape, dtype=float)
    if prepared_custom_solar is not None:
        return prepared_custom_solar.interpolate(wavenumber)
    solar_spectral_irradiance, solar_wavenumber = (
        utilities.load_solar_spectrum_Gueymard20018()
    )
    solar_spectral_irradiance = np.interp(
        wavenumber,
        solar_wavenumber[::-1],
        solar_spectral_irradiance[::-1],
    )
    return utilities.scale_solar_spectrum(solar_spectral_irradiance, year_day)


def _set_spectral_boundary_inputs(
    disort_model,
    spectral_index,
    surface_albedo,
    solar_spectral_irradiance,
    *,
    sun,
):
    """Set per-wavenumber Lambertian albedo and beam-normal FBEAM.

    Args:
        disort_model (DISORT): Configured DISORT wrapper.
        spectral_index (int): Current computational-grid index.
        surface_albedo (float | numpy.ndarray): Scalar fast path or spectral values.
        solar_spectral_irradiance (array-like): Prepared FBEAM spectrum.
        sun (bool): Pass the prepared beam when true and exactly zero otherwise.
    """
    if not isinstance(surface_albedo, float):
        disort_model.set_albedo(surface_albedo[spectral_index])
    disort_model.set_fbeam(
        solar_spectral_irradiance[spectral_index] if sun else 0
    )


def _construct_and_validate_optical_layers(values, computational_wavenumber_cm_inverse):
    """Construct and validate every configured optical-layer object.

    All objects are validated before any one of them calculates optical
    properties. Cross-layer geometry has already been checked by the top-level
    schema, while this object-level pass enforces the same resolved bounds used
    during insertion.

    Args:
        values (Mapping): Validated runner inputs.
        computational_wavenumber_cm_inverse (array-like): Complete model grid.

    Returns:
        tuple[dict, dict, dict]: Mie, prescribed, and grey-body layer mappings.

    Raises:
        ValueError: If resolved layer objects overlap or share a boundary.
    """
    mie_layers = {}
    for layer_name, layer_inputs in (values.get("scat_lyrs_inputs") or {}).items():
        mie_layer = layer.MieLayer()
        mie_layer.set_input_from_dict(dict(layer_inputs))
        mie_layers[layer_name] = mie_layer

    prescribed_layers = {}
    for layer_name, layer_inputs in (values.get("prescribed_lyrs_inputs") or {}).items():
        prescribed_layer = layer.PrescribedOpticalLayer()
        prescribed_layer.set_input_from_dict(dict(layer_inputs))
        prescribed_layers[layer_name] = prescribed_layer

    grey_body_layers = {}
    for layer_name, layer_inputs in (values.get("gbc_lyrs_inputs") or {}).items():
        grey_body_layer = layer.GreyBodyCloud()
        grey_body_layer.set_input_from_dict(dict(layer_inputs))
        grey_body_layers[layer_name] = grey_body_layer

    for mie_layer in mie_layers.values():
        mie_layer.validate_inputs()
    for prescribed_layer in prescribed_layers.values():
        prescribed_layer.validate_inputs(
            computational_wavenumber_cm_inverse,
            scattering_block_size=values.get("scattering_block_size", 10000),
        )
    for grey_body_layer in grey_body_layers.values():
        grey_body_layer.validate_inputs()

    resolved_layers = [
        (name, configured)
        for group in (mie_layers, prescribed_layers, grey_body_layers)
        for name, configured in group.items()
    ]
    resolved_layers.sort(key=lambda item: (item[1].alt_low, item[1].alt_upp, item[0]))
    for first_index, (first_name, first_layer) in enumerate(resolved_layers):
        for second_name, second_layer in resolved_layers[first_index + 1 :]:
            if second_layer.alt_low > first_layer.alt_upp:
                break
            raise ValueError(
                f"Optical layers {first_name!r} and {second_name!r} overlap or "
                "share a vertical boundary."
            )
    return mie_layers, prescribed_layers, grey_body_layers


def _calculate_optical_layers(mie_layers, grey_body_layers, retain_phase_functions):
    """Calculate validated Mie and grey-body optical properties."""
    for mie_layer in mie_layers.values():
        mie_layer.calculate_op(validate=False)
        mie_layer.add_op_calc_output()
        if not retain_phase_functions:
            mie_layer.discard_phase_function()
    for grey_body_layer in grey_body_layers.values():
        grey_body_layer.calculate_op(validate=False)


_GREY_BODY_LAYER_ATTRIBUTES = (
    "name",
    "low_spc",
    "upp_spc",
    "spec_units",
    "res",
    "center_alt",
    "thick",
    "alt_low",
    "alt_upp",
    "emis",
    "inp_tau",
)


def _prepare_grey_body_layers(values):
    """Construct and calculate every configured grey-body cloud layer.

    Args:
        values (Mapping): Validated runner input mapping.

    Returns:
        dict[str, GreyBodyCloud]: Calculated layers keyed by their configured names.
    """
    grey_body_layers = {}
    for layer_name, layer_inputs in (values.get("gbc_lyrs_inputs") or {}).items():
        grey_body_layer = layer.GreyBodyCloud()
        grey_body_layer.set_input_from_dict(dict(layer_inputs))
        grey_body_layer.calculate_op()
        grey_body_layers[layer_name] = grey_body_layer
    return grey_body_layers


def _grey_body_effective_parameters(grey_body_layers):
    """Return serializable effective inputs, including calculated layer bounds."""
    return {
        layer_name: {"layer_type": "grey_body"}
        | {
            attribute: getattr(grey_body_layer, attribute)
            for attribute in _GREY_BODY_LAYER_ATTRIBUTES
            if hasattr(grey_body_layer, attribute)
        }
        for layer_name, grey_body_layer in grey_body_layers.items()
    }


def _interpolate_grey_body_optical_depths(grey_body_layers, wavelengths):
    """Interpolate non-scattering cloud depths onto the actual RFM grid."""
    return {
        layer_name: grey_body_layer.interpolate_optical_depth(wavelengths)
        for layer_name, grey_body_layer in grey_body_layers.items()
    }


def _add_grey_body_optical_depth(
    gas_optical_depth,
    tracked_layers,
    grey_body_optical_depths,
    spectral_index,
):
    """Add grey-body absorption to a copy of one RFM optical-depth row."""
    combined = np.asarray(gas_optical_depth, dtype=float).copy()
    for layer_name, optical_depth in grey_body_optical_depths.items():
        combined[tracked_layers.index(layer_name)] += optical_depth[spectral_index]
    return combined


def _calculate_output_utau(
    output_format,
    output_values,
    altitude_rows,
    layer_optical_depths,
    first_retained_layer,
):
    """Convert configured output levels to the optical depths used by DISORT.

    Args:
        output_format (str): ``"altitude"`` or direct optical-depth ``"tau"``.
        output_values (array-like): Resolved altitude or optical-depth outputs.
        altitude_rows (sequence[int | None] | None): Atmospheric rows matching
            altitude outputs, with ``None`` representing top of atmosphere.
        layer_optical_depths (array-like): Top-to-bottom total layer depths.
        first_retained_layer (int): Index of the first non-truncated layer.

    Returns:
        list[float]: Optical depths measured from the retained atmospheric top.

    Raises:
        ValueError: If a requested output lies below the retained atmosphere.
    """
    layer_optical_depths = np.asarray(layer_optical_depths, dtype=float)
    retained = layer_optical_depths[first_retained_layer:]

    if output_format == "tau":
        utau = np.asarray(output_values, dtype=float).copy()
    else:
        cumulative = np.cumsum(layer_optical_depths)
        discarded = float(layer_optical_depths[:first_retained_layer].sum())
        utau = np.asarray(
            [
                0.0 if row is None else max(float(cumulative[row]) - discarded, 0.0)
                for row in altitude_rows
            ]
        )

    total_optical_depth = float(retained.sum())
    tolerance = max(1e-12, abs(total_optical_depth) * 1e-10)
    if np.any(utau > total_optical_depth + tolerance):
        raise ValueError(
            "Requested output optical depth exceeds the atmospheric optical depth "
            f"({total_optical_depth:g}) at this wavenumber."
        )
    return np.minimum(utau, total_optical_depth).tolist()


def _write_spectral_text(filename, wavenumbers, values, value_name):
    """Write every spectrum in a four-dimensional SRFM output array.

    Args:
        filename (path-like): Destination text file.
        wavenumbers (array-like): Spectral coordinate in cm-1.
        values (array-like): Values ordered as spectral, polar, level, azimuth.
        value_name (str): Descriptive column heading.

    Returns:
        None: The formatted spectra are written to disk.

    Raises:
        ValueError: If the value array does not have the expected shape.
    """
    values = np.asarray(values)
    if values.ndim != 4 or values.shape[0] != len(wavenumbers):
        raise ValueError("Text output requires an array shaped (spectral, polar, level, azimuth).")

    _, num_polar, num_levels, num_azimuthal = values.shape
    spectra_per_wavenumber = num_polar * num_levels * num_azimuthal
    polar, level, azimuthal = np.indices(
        (num_polar, num_levels, num_azimuthal)
    )
    indices = (polar.ravel(), level.ravel(), azimuthal.ravel())

    with open(filename, "w", encoding="utf-8") as output_file:
        output_file.write(
            "# Wavenumber (cm-1), polar angle index, output level index, "
            f"azimuthal angle index, {value_name}\n"
        )
        for start in range(0, values.shape[0], 10_000):
            end = min(start + 10_000, values.shape[0])
            block_size = end - start
            block = np.column_stack(
                (
                    np.repeat(wavenumbers[start:end], spectra_per_wavenumber),
                    np.tile(indices[0], block_size),
                    np.tile(indices[1], block_size),
                    np.tile(indices[2], block_size),
                    values[start:end].reshape(-1),
                )
            )
            np.savetxt(output_file, block, fmt=["%.4f", "%d", "%d", "%d", "%.8e"])


def _create_netcdf_spectral_dimensions(nc_file, model):
    """Create dimensions and coordinates shared by retained SRFM outputs.

    Args:
        nc_file (netCDF4.Dataset): Open output dataset.
        model (SRFM): Model containing geometry and spectral coordinates.

    Returns:
        tuple[str, str, str, str]: NetCDF dimension names in native output order.
    """
    num_wavenumbers = len(model.wvnm)
    num_polar = len(model.output_polar_angles)
    num_levels = len(model.output_values)
    num_azimuthal = len(model.output_azimuthal_angles)
    nc_file.createDimension("wavenumber", num_wavenumbers)
    nc_file.createDimension("output_polar_angle", num_polar)
    nc_file.createDimension("output_level", num_levels)
    nc_file.createDimension("output_azimuthal_angle", num_azimuthal)

    wavenumber = nc_file.createVariable("wavenumber", "f8", ("wavenumber",))
    wavenumber.units = "cm-1"
    wavenumber.long_name = "Wavenumber"
    wavenumber[:] = model.wvnm

    output_level = nc_file.createVariable("output_level", "f8", ("output_level",))
    output_level[:] = model.output_values
    if model.output_format == "altitude":
        output_level.units = "km"
        output_level.long_name = "Matched output altitude"
    else:
        output_level.units = "1"
        output_level.long_name = "Output optical depth"

    polar_index = nc_file.createVariable(
        "output_polar_angle_index", "i4", ("output_polar_angle",)
    )
    polar_index.long_name = "Index of output polar angle"
    polar_index[:] = np.arange(num_polar)

    azimuthal_index = nc_file.createVariable(
        "output_azimuthal_angle_index", "i4", ("output_azimuthal_angle",)
    )
    azimuthal_index.long_name = "Index of output azimuthal angle"
    azimuthal_index[:] = np.arange(num_azimuthal)

    return (
        "wavenumber",
        "output_polar_angle",
        "output_level",
        "output_azimuthal_angle",
    )


def _write_retained_netcdf_outputs(nc_file, model, retained_outputs):
    """Write explicitly retained SRFM outputs to an open NetCDF dataset.

    The SRFM arrays use three native DISORT layouts: output-level fields,
    output-polar-angle fields, and the full user-angle radiance field.  The
    shared NetCDF dimensions must already have been created by
    :func:`_create_netcdf_spectral_dimensions`.

    Args:
        nc_file (netCDF4.Dataset): Open output dataset.
        model (SRFM): Model containing the retained arrays.
        retained_outputs (collection[str] | None): Values explicitly selected
            by the user. ``"rad"`` and ``"radiance"`` are treated as aliases
            for ``"uu"``.

    Returns:
        None: Selected arrays are written as compressed NetCDF variables.

    Raises:
        RuntimeError: If a selected array is absent from the model.
        ValueError: If a selected name or array shape is unsupported.
    """
    if retained_outputs is None:
        return

    spectrum_dimensions = (
        "wavenumber",
        "output_polar_angle",
        "output_level",
        "output_azimuthal_angle",
    )
    level_dimensions = ("wavenumber", "output_level")
    polar_dimensions = ("wavenumber", "output_polar_angle")
    output_metadata = {
        "rfldir": (level_dimensions, "Direct-beam downward flux"),
        "rfldn": (level_dimensions, "Diffuse downward flux"),
        "flup": (level_dimensions, "Diffuse upward flux"),
        "dfdt": (level_dimensions, "Flux-divergence derivative"),
        "uavg": (level_dimensions, "Mean intensity"),
        "uu": (spectrum_dimensions, "User-angle radiance"),
        "albmed": (polar_dimensions, "Medium albedo"),
        "trnmed": (polar_dimensions, "Medium transmissivity"),
        "bbt": (spectrum_dimensions, "Brightness temperature"),
    }

    aliases = {"rad": "uu", "radiance": "uu"}
    selected = {aliases.get(name, name) for name in retained_outputs}
    unknown = selected - output_metadata.keys()
    if unknown:
        raise ValueError(
            "Unknown retained NetCDF output name(s): "
            + ", ".join(sorted(unknown))
        )

    for output_name in sorted(selected):
        if not hasattr(model, output_name):
            raise RuntimeError(
                f"Retained output {output_name!r} is not available on the SRFM model."
            )

        dimensions, long_name = output_metadata[output_name]
        values = np.asarray(getattr(model, output_name))
        expected_shape = tuple(len(nc_file.dimensions[name]) for name in dimensions)
        if values.shape != expected_shape:
            raise ValueError(
                f"Retained output {output_name!r} has shape {values.shape}, "
                f"but dimensions {dimensions} require {expected_shape}."
            )

        variable = nc_file.createVariable(
            output_name, "f8", dimensions, zlib=True, complevel=4
        )
        variable.long_name = long_name
        if output_name == "bbt":
            variable.units = "K"
        elif output_name == "uu":
            variable.units = "W m-2 sr-1 cm"
        variable[:] = values


def _compact_metadata(value):
    """Return JSON-safe provenance without copying or embedding large arrays."""
    if isinstance(value, Mapping):
        compact = {}
        for key, item in value.items():
            if key in {"grid", "values"} and isinstance(
                item, (np.ndarray, list, tuple)
            ):
                array = np.asarray(item)
                compact[key] = {
                    "shape": list(array.shape),
                    "dtype": str(array.dtype),
                    "stored_as": "NetCDF variable",
                }
            else:
                compact[key] = _compact_metadata(item)
        return compact
    if isinstance(value, np.ndarray):
        return {
            "shape": list(value.shape),
            "dtype": str(value.dtype),
            "stored_as": "NetCDF variable",
        }
    if isinstance(value, (list, tuple)):
        if len(value) > 100:
            array = np.asarray(value)
            return {
                "shape": list(array.shape),
                "dtype": str(array.dtype),
                "stored_as": "NetCDF variable",
            }
        return [_compact_metadata(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    return value


def _write_srfm_netcdf(
    filename,
    model,
    retained_outputs,
    effective_params,
    computational_wavenumber_cm_inverse,
    particle_layers,
    grey_body_layers=None,
    surface_albedo=None,
    solar_spectral_irradiance=None,
    nmom=None,
    scattering_block_size=10000,
):
    """Write one complete SRFM NetCDF file.

    The file contains every explicitly retained SRFM result, shared output
    coordinates, and the existing scattering-layer optical-depth metadata.

    Args:
        filename (path-like): Destination NetCDF filename.
        model (SRFM): Completed model containing retained result arrays.
        retained_outputs (collection[str]): User-selected result names.
        effective_params (Mapping): Effective validated run configuration.
        computational_wavenumber_cm_inverse (array-like): Computational grid in
            ``cm-1``.
        particle_layers (Mapping[str, Layer]): Mie and prescribed particle layers.
        grey_body_layers (Mapping[str, GreyBodyCloud] | None): Absorbing clouds.
        surface_albedo (float | array-like | None): Albedo used on the
            computational grid.
        solar_spectral_irradiance (array-like | None): Beam-normal FBEAM spectrum
            used on the computational grid.
        nmom (int | None): Highest particle-moment order to serialize.
        scattering_block_size (int): Maximum interpolation block size.

    Returns:
        None: The NetCDF file is written and closed.
    """
    computational_wavenumber_cm_inverse = np.asarray(
        computational_wavenumber_cm_inverse, dtype=float
    )
    particle_layers = particle_layers or {}
    grey_body_layers = grey_body_layers or {}
    layer_names = list(particle_layers) + list(grey_body_layers)
    serialized_source = dict(effective_params)
    if isinstance(serialized_source.get("driver_inputs"), Mapping):
        driver_inputs = dict(serialized_source["driver_inputs"])
        if "spectral" in driver_inputs:
            driver_inputs["spectral"] = str(driver_inputs["spectral"])
        serialized_source["driver_inputs"] = driver_inputs
    serialized_params = _compact_metadata(serialized_source)

    with Dataset(filename, "w", format="NETCDF4") as nc_file:
        nc_file.description = "SRFM output."
        nc_file.history = (
            f"Created {datetime.datetime.now().strftime('%Y-%m-%d')}"
        )
        nc_file.srfm_params = json.dumps(
            serialized_params,
            default=utilities.json_handler,
        )

        nc_file.createDimension(
            "wavenumber_op", computational_wavenumber_cm_inverse.size
        )
        nc_file.createDimension("layer", len(layer_names))
        _create_netcdf_spectral_dimensions(nc_file, model)
        _write_retained_netcdf_outputs(nc_file, model, retained_outputs)

        optical_wavenumber = nc_file.createVariable(
            "wavenumber_op", "f8", ("wavenumber_op",)
        )
        optical_wavenumber.units = "cm-1"
        optical_wavenumber.long_name = "Computational optical-property wavenumber"
        optical_wavenumber[:] = computational_wavenumber_cm_inverse

        layer_variable = nc_file.createVariable("layer_names", str, ("layer",))
        layer_variable[:] = np.asarray(layer_names, dtype=object)
        layer_type = nc_file.createVariable("layer_type", str, ("layer",))
        layer_type.long_name = "Configured optical-layer representation"
        layer_type[:] = np.asarray(
            [
                (
                    "mie"
                    if isinstance(particle_layers[name], layer.MieLayer)
                    else "prescribed"
                )
                if name in particle_layers
                else "grey_body"
                for name in layer_names
            ],
            dtype=object,
        )

        tau = nc_file.createVariable(
            "tau",
            "f8",
            ("layer", "wavenumber_op"),
            zlib=True,
            complevel=4,
        )
        tau.long_name = "Optical depth per layer and wavenumber"

        if nmom is None:
            nmom = max(
                (configured.required_nmom for configured in particle_layers.values()),
                default=0,
            )
        if particle_layers:
            nc_file.createDimension("moment", nmom + 1)
            moment_order = nc_file.createVariable("moment", "i4", ("moment",))
            moment_order.long_name = "Legendre moment order"
            moment_order[:] = np.arange(nmom + 1)
            single_scattering_albedo = nc_file.createVariable(
                "single_scattering_albedo",
                "f8",
                ("layer", "wavenumber_op"),
                zlib=True,
                complevel=4,
            )
            single_scattering_albedo.long_name = "Particle single-scattering albedo"
            normalized_moments = nc_file.createVariable(
                "normalized_legendre_moment",
                "f8",
                ("layer", "wavenumber_op", "moment"),
                zlib=True,
                complevel=4,
            )
            normalized_moments.long_name = "Normalized particle phase-function moment"

        for layer_index, layer_name in enumerate(layer_names):
            if layer_name in particle_layers:
                configured_layer = particle_layers[layer_name]
                for start in range(
                    0,
                    computational_wavenumber_cm_inverse.size,
                    scattering_block_size,
                ):
                    stop = min(
                        start + scattering_block_size,
                        computational_wavenumber_cm_inverse.size,
                    )
                    properties = configured_layer.column_optical_properties(
                        computational_wavenumber_cm_inverse[start:stop], nmom=nmom
                    )
                    tau[layer_index, start:stop] = properties[
                        "particle_optical_depth"
                    ]
                    single_scattering_albedo[layer_index, start:stop] = properties[
                        "single_scattering_albedo"
                    ]
                    normalized_moments[layer_index, start:stop, :] = properties[
                        "normalized_legendre_coefficients"
                    ]
            else:
                wavelengths = 1.0e4 / computational_wavenumber_cm_inverse
                tau[layer_index, :] = grey_body_layers[
                    layer_name
                ].interpolate_optical_depth(wavelengths)
                if particle_layers:
                    single_scattering_albedo[layer_index, :] = 0.0
                    normalized_moments[layer_index, :, :] = 0.0
                    normalized_moments[layer_index, :, 0] = 1.0

        if surface_albedo is not None:
            albedo_values = np.asarray(surface_albedo, dtype=float)
            if albedo_values.ndim == 0:
                albedo_values = np.full(
                    computational_wavenumber_cm_inverse.shape, albedo_values.item()
                )
            albedo_variable = nc_file.createVariable(
                "surface_albedo", "f8", ("wavenumber_op",), zlib=True, complevel=4
            )
            albedo_variable.long_name = "Lambertian surface albedo"
            albedo_variable.units = "1"
            albedo_variable[:] = albedo_values

        if solar_spectral_irradiance is not None:
            solar_variable = nc_file.createVariable(
                "solar_spectral_irradiance",
                "f8",
                ("wavenumber_op",),
                zlib=True,
                complevel=4,
            )
            solar_variable.long_name = "DISORT beam-normal FBEAM spectral irradiance"
            solar_variable.units = "W m-2 (cm-1)-1"
            solar_variable[:] = np.asarray(solar_spectral_irradiance, dtype=float)


def _plot_spectral_outputs(
    model, results_folder, output_name, show_plots, filename_prefix="base_plot"
):
    """Plot every polar-angle, output-level, and azimuthal-angle spectrum."""
    if output_name == "bbt":
        values = model.bbt
        y_label = "Brightness temperature (K)"
    elif output_name == "rad":
        values = model.uu
        y_label = r"Radiance (W m$^{-2}$ sr$^{-1}$ cm)"
    else:
        raise ValueError("Plot type not recognized. Please use 'bbt' or 'rad'.")

    level_name = (
        "altitude (km)" if model.output_format == "altitude" else "optical depth"
    )
    for polar in range(values.shape[1]):
        for level in range(values.shape[2]):
            for azimuthal in range(values.shape[3]):
                plt.figure(figsize=(11.7, 8.4))
                plt.rcParams.update({"font.size": 12})
                plt.plot(
                    model.wvnm,
                    values[:, polar, level, azimuthal],
                    label="SRFM",
                    color="tab:blue",
                )
                plt.xlabel(r"Wavenumbers (cm$^{-1}$)")
                plt.ylabel(y_label)
                plt.title(
                    f"Output polar angle: {model.output_polar_angles[polar]:g}°\n"
                    f"Output {level_name}: {model.output_values[level]:g}\n"
                    f"Output azimuthal angle: {model.output_azimuthal_angles[azimuthal]:g}°"
                )
                plt.legend()
                filename = (
                    f"{filename_prefix}_{polar}_{level}_{azimuthal}.png"
                )
                plt.savefig(os.path.join(results_folder, filename))
                if show_plots:
                    plt.ion()
                    plt.show()
                else:
                    plt.close()


def _plot_retained_spectral_outputs(
    model, results_folder, retained_outputs, show_plots
):
    """Create base plots for every retained primary spectral output.

    When both brightness temperature and radiance are retained, the output
    name is included in each filename so the two plot sets cannot overwrite
    one another.
    """
    selected = {
        "uu" if name in {"rad", "radiance"} else name
        for name in retained_outputs
    }
    plot_outputs = [
        output_name
        for output_name, retained_name in (("bbt", "bbt"), ("rad", "uu"))
        if retained_name in selected
    ]
    if not plot_outputs:
        warnings.warn(
            "base_plots is True, but retain_outputs contains neither bbt nor "
            "radiance; no base plots were created.",
            stacklevel=2,
        )
        return

    for output_name in plot_outputs:
        filename_prefix = (
            "base_plot"
            if len(plot_outputs) == 1
            else f"base_plot_{output_name}"
        )
        _plot_spectral_outputs(
            model,
            results_folder,
            output_name,
            show_plots,
            filename_prefix=filename_prefix,
        )


@utilities.show_runtime
def run_srfm(inp):
    """Run the generic scattering reference forward model.

    RFM optical depths are consumed in a compact spectral-first array and particle
    properties are interpolated in bounded blocks before individual DISORT calls.

    Generic srfm  run.

    Args:
        inp (obj): Instance of inputs.Inputs.

    Returns:
        model_SRFM (obj): Instance of forward_model.SRFM.

    Raises:
        TypeError: If ``inp`` does not provide a ``values`` mapping.
        InputValidationError: If configured inputs are invalid.
        RuntimeError: If RFM capture or a native model calculation fails.
        ValueError: If spectral, atmospheric, or output geometry is inconsistent.

    """
    if not hasattr(inp, "values"):
        raise TypeError("run_srfm expects an Inputs-like object with a values mapping.")
    inp.values = validate_srfm_inputs(inp.values)

    ########################################################################################
    # set final grid to interpolate to
    ########################################################################################
    # this seems to be the most robust way of generating a grid (both np.arange and linspace are prone to failing)
    npts = (
        int(
            np.floor(
                (inp.values["fin_wvnmhi"] - inp.values["fin_wvnmlo"])
                / inp.values["fin_res"]
            )
        )
        + 1
    )  # expected number of points in the grid
    fin_grid = inp.values["fin_wvnmlo"] + np.arange(npts) * inp.values["fin_res"]

    if "date" in inp.values:
        if isinstance(inp.values["date"], datetime.datetime):
            date = inp.values["date"]
        elif isinstance(inp.values["date"], tuple) and len(inp.values["date"]) == 3:
            date = datetime.datetime(*inp.values["date"])
        else:
            raise ValueError(f"""date must be a tuple of len 3 or a 
                datetime.datetime object. The current date is 
                {inp.values["date"]} and is a type 
                {type(inp.values["date"])}.""")

    else:
        date = datetime.datetime(
            2025, 3, 23
        )  # default date, approx spring equinox, avg Earth-Sun dist
    year_day = date.timetuple().tm_yday

    ########################################################################################
    # specify spectral calculation grid
    ########################################################################################
    spec_res = inp.values["spc_res"]  # model spectral resolution,[spec_units]
    low_spc = inp.values["spc_wvnmlo"]  # model start wavenumber (lower), [spec_units]
    upp_spc = inp.values["spc_wvnmhi"]  # model end wavenumber (upper), [spec_units]
    spec_units = inp.values["spc_units"]  # accepted values "cm-1", "um", "nm"

    RFM_wvnm, wvls = utilities.calc_grids(low_spc, upp_spc, spec_res, spec_units)

    # Load and validate every new spectral source and every configured layer before
    # creating output files or invoking an optical/native model calculation.
    prepared_albedo, prepared_custom_solar = _prepare_boundary_spectral_fields(
        inp.values, RFM_wvnm
    )
    scat_lyrs, prescribed_lyrs, gbc_lyrs = _construct_and_validate_optical_layers(
        inp.values, RFM_wvnm
    )
    _calculate_optical_layers(
        scat_lyrs,
        gbc_lyrs,
        inp.values.get("retain_phase_functions", False),
    )
    particle_lyrs = {**scat_lyrs, **prescribed_lyrs}

    # Keep all run-generated files outside the installed package tree.
    os.makedirs(inp.values["results_fldr"], exist_ok=True)
    rfm_grid_fname = rfm_functions.construct_rfm_grid_file(
        RFM_wvnm, filename="grid.spc", rfm_fldr=inp.values["results_fldr"]
    )

    # Prepare dict with layer parameters to be saved in the output
    layer_attrs = (
        "name",
        "low_spc",
        "upp_spc",
        "spec_units",
        "res",
        "mass_loading",
        "n",
        "r",
        "s",
        "rho",
        "s_a_den",
        "v_den",
        "dist_type",
        "comp",
        "center_alt",
        "thick",
        "alt_low",
        "alt_upp",
        "radii",
        "eta",
        "phase_quad_N",
        "phase_quad_type",
        "radii_quad_type",
        "leg_coeffs",
        "leg_coeffs_type",
        "multiprocess",
    )

    effective_params = dict(inp.values)
    effective_params["scat_lyrs_inputs"] = {
        lyr: {"layer_type": "mie"}
        | {
            attr: getattr(scat_lyrs[lyr], attr)
            for attr in layer_attrs
            if hasattr(scat_lyrs[lyr], attr)
        }
        for lyr in scat_lyrs
    }
    effective_params["prescribed_lyrs_inputs"] = {
        lyr: prescribed_lyrs[lyr].source_metadata() for lyr in prescribed_lyrs
    }
    effective_params["gbc_lyrs_inputs"] = _grey_body_effective_parameters(gbc_lyrs)
    effective_params["albedo"] = (
        prepared_albedo
        if isinstance(prepared_albedo, float)
        else prepared_albedo.provenance()
    )
    if prepared_custom_solar is not None:
        effective_params["solar_spectrum"] = prepared_custom_solar.provenance()

    ########################################################################################
    # prepare atmospheric layer structure
    ########################################################################################

    # define some requested output levels
    if "levels" in inp.values.keys():
        levels = inp.values["levels"]
    else:
        levels = [
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            8.0,
            9.0,
            10.0,
            11.0,
            12.0,
            13.0,
            14.0,
            15.0,
            16.0,
            17.0,
            18.0,
            19.0,
            20.0,
            21.0,
            22.0,
            23.0,
            24.0,
            25.0,
            27.5,
            30.0,
            32.5,
            35.0,
            37.5,
            40.0,
            42.5,
            45.0,
            47.5,
            50.0,
            55.0,
            60.0,
            65.0,
            70.0,
            75.0,
            80.0,
            85.0,
            90.0,
            95.0,
            100.0,
            115.0,
            120.0,
        ]

    # define tracker array for atmospheric structure
    track_lev = [None for i in levels]

    # add upper and lower particle layer boundaries, delete any levels "within" the layer
    for lyr in scat_lyrs:
        levels, track_lev = utilities.add_lyr_from_Layer(
            lev=levels, track_lev=track_lev, new_lyr=scat_lyrs[lyr]
        )

    for lyr in prescribed_lyrs:
        levels, track_lev = utilities.add_lyr_from_Layer(
            lev=levels, track_lev=track_lev, new_lyr=prescribed_lyrs[lyr]
        )

    for lyr in gbc_lyrs:
        levels, track_lev = utilities.add_lyr_from_Layer(
            lev=levels, track_lev=track_lev, new_lyr=gbc_lyrs[lyr]
        )

    # convert the tracking levels array to a tracking layers array
    track_lyr = utilities.track_lev_to_track_lyr(track_lev)
    track_lyr = track_lyr[::-1]
    
#    # insert "aerosol" into RFM: prepare xsc file and atm file
#    rfm_prf = rfm_functions.read_atm_file(inp.values["driver_inputs"]["atmosphere"][1])
#    hgt_key = [i for i in rfm_prf if i.lower().startswith("hgt ")] 
#    tem_key = [i for i in rfm_prf if i.lower().startswith("pre ")]
#    pre_key = [i for i in rfm_prf if i.lower().startswith("tem ")]
#    
#    xscs = {} # dictionary with extinxtions for xsc file, keys are just numbers
#    
#    for i,lyr in enumerate(scat_lyrs):
#        xscs[f"{0}l"] = {}
#        xscs[f"{0}l"]["molec"] = "Aerosol"
#        xscs[f"{0}l"][tem_key] = np.interp(lyr.alt_low,rfm_prf[hgt_key],rfm_prf[tem_key])
#        xscs[f"{0}l"][pre_key] = np.exp(np.interp(lyr.alt_low,rfm_prf[hgt_key], np.log(np.maximum(rfm_prf[pre_key],tiny))))
#        xscs[f"{0}l"]["low_spc"] = lyr.low_spc
#        xscs[f"{0}l"]["upp_spc"] = lyr.upp_spc 
#        xscs[f"{0}l"]["npts"] = (int(np.floor((lyr.upp_spc - lyr.low_spc)/ lyr.res))+ 1)  # expected number of points in the grid
#        xscs[f"{0}l"]["xsc"] = lyr.beta_ext # units m-1
#        
#        xscs[f"{0}u"] = {}
#        xscs[f"{0}u"]["molec"] = "Aerosol"
#        xscs[f"{0}u"][tem_key] = np.interp(lyr.alt_upp,rfm_prf[hgt_key],rfm_prf[tem_key])
#        xscs[f"{0}u"][pre_key] = np.exp(np.interp(lyr.alt_upp,rfm_prf[hgt_key], np.log(np.maximum(rfm_prf[pre_key],tiny))))
#        xscs[f"{0}u"]["low_spc"] = lyr.low_spc
#        xscs[f"{0}u"]["upp_spc"] = lyr.upp_spc 
#        xscs[f"{0}u"]["npts"] = (int(np.floor((lyr.upp_spc - lyr.low_spc)/ lyr.res))+ 1)  # expected number of points in the grid
#        xscs[f"{0}u"]["xsc"] = lyr.beta_ext # units m-1
#    
#    # write xsc file
#    rfm_functions.write_xsc_file(xscs, filename=os.path.join(inp.values["results_fldr"],"aerosol.xsc"))
#    
#    aslprf = {}
#    aslprf["HGT (km)"] = levels
#    aslmask = [1 if track_lev[levels.index(i)] != None else 0 for i in levels]
#    rfm_functions.write_atm_file(aslmask,filename=f"{inp.values['results_fldr']}/aerosol.atm")

    ########################################################################################
    # prepare and call RFM
    ########################################################################################
    # precalculate angles
    zen_rad = np.deg2rad(inp.values["zen"])  # zenith angle to rad
    zen_cos = np.cos(zen_rad).item()  # cosine of zenith angle
    zen_sec = 1 / zen_cos  # secant of zenith angle

    # RFM global config
    rfm_config = inp.values["rfm_config"]

    # RFM driver table
    driver_inputs = copy.deepcopy(inp.values["driver_inputs"])
    driver_inputs["spectral"] = (rfm_helper.SpectralFile(rfm_grid_fname),)
    driver_inputs["lev"] = tuple(str(val) for val in levels)
    driver_inputs["tangent"] = (str(zen_sec),)
    
    if "btemp" in inp.values:
        driver_inputs["sfc"] = (f"TEMSFC={inp.values['btemp']}",)


#    driver_inputs["atmosphere"] = list(driver_inputs["atmosphere"])
#    driver_inputs["atmosphere"].append(f"{inp.values['results_fldr']}/aerosol.atm")
#    driver_inputs["atmosphere"] = tuple(driver_inputs["atmosphere"])

#    driver_inputs["xsc"] = list(driver_inputs["xsc"])
#    driver_inputs["xsc"].append(f"{inp.values['results_fldr']}/aerosol.xsc")
#    driver_inputs["xsc"] = tuple(driver_inputs["xsc"])

    # initialize RFM model class
    model_RFM = forward_model.RFM()

    # print current status:
    print(model_RFM.status)

    # run rfm
    rfm_run_result = rfm_helper.rfm_main(
        configuration=rfm_config,
        driver_inputs=driver_inputs,
        levels=levels,
        rfm_out_fldr=inp.values["results_fldr"],
    )
    # Keep only compact numerical optical-depth data on the main path.
    model_RFM.rfm_run_result = rfm_run_result
    if rfm_run_result.optical_depth_grid is None:
        raise RuntimeError("RFM capture did not return a compact optical-depth grid.")
    model_RFM.rfm_output = rfm_run_result.optical_depth_grid

    # print current status:
    model_RFM.status = "RFM completed"
    print(model_RFM.status)

    RFM_wvnm = model_RFM.rfm_output.wavenumber
    wvls = (1.0 / RFM_wvnm) * 1e4
    if isinstance(prepared_albedo, float):
        surface_albedo = prepared_albedo
    else:
        surface_albedo = prepared_albedo.interpolate(RFM_wvnm)
    grey_body_optical_depths = _interpolate_grey_body_optical_depths(
        gbc_lyrs, wvls
    )

    ########################################################################################
    # prepare DISORT common variables
    ########################################################################################

    # initialize DISORT model class
    model_DISORT = forward_model.DISORT(retain_history=False)

    # check if number of columns and wavelengths match
    if model_RFM.rfm_output.differential_tau.shape[0] != len(wvls):
        raise ValueError(
            f"Number of RFM and scattering wavelengths don't match "
            f"({model_RFM.rfm_output.differential_tau.shape[0]} optical-depth rows "
            f"vs {len(wvls)} scattering wavelengths)."
        )
    output_values, altitude_rows = _resolve_output_geometry(
        inp.values["out_fmt"],
        inp.values["out"],
        inp.values["out_toa"],
        model_RFM.rfm_output,
    )

    # set disort_input parameters common to all loop iterations
    # these need to be set first:
    nmom = inp.values["nmom"]
    for particle_layer in particle_lyrs.values():
        nmom = max(nmom, particle_layer.required_nmom)

    model_DISORT.set_maxcmu(inp.values["maxcmu"])

    model_DISORT.set_maxmom(nmom)
    if nmom < model_DISORT.disort_input["maxcmu"]:
        model_DISORT.set_maxmom(model_DISORT.disort_input["maxcmu"])
    nmom = model_DISORT.disort_input["maxmom"]

    model_DISORT.set_maxumu(inp.values["maxumu"])
    model_DISORT.set_maxphi(inp.values["maxphi"])
    # maxulv is an internal DISORT dimension derived from the resolved output
    # geometry rather than a public driver-table input.
    model_DISORT.set_maxulv(len(output_values))
    effective_params["maxulv"] = len(output_values)

    # now the rest
    model_DISORT.set_usrang(inp.values["usrang"])
    model_DISORT.set_usrtau(inp.values["usrtau"])
    model_DISORT.set_ibcnd(inp.values["ibcnd"])
    model_DISORT.set_onlyfl(inp.values["onlyfl"])
    # model_DISORT.set_prnt([True, True, True, False, False])
    model_DISORT.set_prnt(inp.values["prnt"])
    model_DISORT.set_plank(inp.values["planck"])
    model_DISORT.set_lamber(inp.values["lamber"])
    model_DISORT.set_deltamplus(inp.values["deltamplus"])
    model_DISORT.set_do_pseudo_sphere(inp.values["do_pseudo_sphere"])

    model_DISORT.set_fisot(inp.values["fisot"])
    if isinstance(surface_albedo, float):
        model_DISORT.set_albedo(surface_albedo)

    model_DISORT.set_temis(inp.values["temis"])
    model_DISORT.set_earth_radius(inp.values.get("earth_radius", 6371.0))
    model_DISORT.set_rhoq(
        np.zeros(
            shape=(
                int(model_DISORT.disort_input["maxcmu"] / 2),
                int(model_DISORT.disort_input["maxcmu"] / 2 + 1),
                int(model_DISORT.disort_input["maxcmu"]),
            )
        )
    )
    model_DISORT.set_rhou(
        np.zeros(
            shape=(
                model_DISORT.disort_input["maxumu"],
                int(model_DISORT.disort_input["maxcmu"] / 2 + 1),
                model_DISORT.disort_input["maxcmu"],
            )
        )
    )
    model_DISORT.set_rho_accurate(
        np.zeros(
            shape=(
                model_DISORT.disort_input["maxumu"],
                model_DISORT.disort_input["maxphi"],
            )
        )
    )
    model_DISORT.set_bemst(
        np.zeros(shape=(int(model_DISORT.disort_input["maxcmu"] / 2)))
    )
    model_DISORT.set_emust(np.zeros(shape=(model_DISORT.disort_input["maxumu"])))
    model_DISORT.set_accur(0)
    ########################################################################################
    # prepare solar spectrum
    ########################################################################################
    solar_spc = _prepare_solar_spectral_irradiance(
        prepared_custom_solar,
        RFM_wvnm,
        year_day,
        sun=inp.values["sun"],
    )
    if inp.values["sun"] == True:
        # get incoming solar beam polar angle for DISORT
        solar_zen_deg = inp.values["sza"]  # solar zenith angle [degrees]
        solar_zen_rad = np.deg2rad(solar_zen_deg)  # solar zenith angle [rad]
        solar_zen_cos = np.cos(
            solar_zen_rad
        ).item()  # cosine of the solar zenith angle, UMU0

        model_DISORT.set_umu0(solar_zen_cos)
        model_DISORT.set_phi0(inp.values["saa"])
    else:
        model_DISORT.set_umu0(1)
        model_DISORT.set_phi0(0)

    # Angles
    model_DISORT.set_phi([inp.values["azi"]])
    model_DISORT.set_umu([zen_cos])

    # initialize disort input arrays for output variables from a single run
    model_DISORT.initialize_disort_output_arrays()

    ########################################################################################
    # set DISORT DISOBRDF variables
    ########################################################################################

    # TBD

    ########################################################################################
    # prepare SRFM common variables
    ########################################################################################
    model_SRFM = forward_model.SRFM()
    model_SRFM.set_wvnm(RFM_wvnm)
    model_SRFM.set_wvls(wvls)
    requested_outputs, runtime_outputs = _resolve_retained_outputs(inp.values)
    model_SRFM.initialize_srfm_output_arrays_from_disort(
        model_DISORT, retain_outputs=runtime_outputs
    )
    model_SRFM.output_format = inp.values["out_fmt"]
    model_SRFM.output_values = output_values
    model_SRFM.output_polar_angles = np.asarray([inp.values["zen"]], dtype=float)
    model_SRFM.output_azimuthal_angles = np.asarray([inp.values["azi"]], dtype=float)

    # track progress
    pct = [
        1,
        2,
        5,
        10,
        20,
        30,
        40,
        50,
        60,
        70,
        80,
        90,
        100,
    ]  # percent completed to report
    pct_val = [RFM_wvnm[int(i * len(RFM_wvnm) / 100 - 1)] for i in pct]

    ########################################################################################
    # set dynamic variables and run DISORT
    ########################################################################################
    scattering_block_size = inp.values.get("scattering_block_size", 10000)
    particle_block = {}
    particle_block_start = 0
    skipped_scattering = {
        layer_name: {"count": 0, "first": None, "last": None}
        for layer_name in particle_lyrs
    }
    for wvl_idx, (wvnm, wvl, tau_g) in enumerate(
        zip(RFM_wvnm, wvls, model_RFM.rfm_output.differential_tau)
    ):
        if particle_lyrs and wvl_idx % scattering_block_size == 0:
            particle_block_start = wvl_idx
            block_stop = min(wvl_idx + scattering_block_size, len(wvls))
            particle_block = _interpolate_particle_block(
                particle_lyrs,
                RFM_wvnm[wvl_idx:block_stop],
                nmom,
            )
        block_index = wvl_idx - particle_block_start

        # track progress
        if wvnm in pct_val:
            val_idx = pct_val.index(wvnm)
            print(f"Running main DISORT loop. {pct[val_idx]}% done...")

        #        model_DISORT.set_header(f"Now starting calculation for {col} cm-1.")
        model_DISORT.set_header(
            inp.values.get("header", "NO HEADER")
        )  # the string "NO HEADER" will cause no printout
        model_DISORT.set_wvnm(wvnm)
        model_DISORT.set_wvl(wvl)

        tau_g = _add_grey_body_optical_depth(
            tau_g,
            track_lyr,
            grey_body_optical_depths,
            wvl_idx,
        )

        # layer optical depths from Rayleigh scattering
        #    tau_R = np.zeros(shape=(tau_g.shape))
        tau_R = utilities.calc_Rayleigh_opt_depths(
            ps=model_RFM.rfm_output.pressure_lower[-1],
            pu=model_RFM.rfm_output.pressure_upper,
            pl=model_RFM.rfm_output.pressure_lower,
            l=wvnm,
        )

        # particle layer optical depths (from particle scattering)
        tau_p = np.zeros(shape=(len(tau_g)))
        for lyr in particle_lyrs:
            tau_p[track_lyr.index(lyr)] = particle_block[lyr][
                "particle_optical_depth"
            ][block_index]

        # particle layer single scatter albedo
        w_p = np.zeros(shape=(len(tau_g)))
        for lyr in particle_lyrs:
            w_p[track_lyr.index(lyr)] = particle_block[lyr][
                "single_scattering_albedo"
            ][block_index]

        dtauc_tot = utilities.calc_tot_dtauc(tau_g=tau_g, tau_R=tau_R, tau_p=tau_p)
        dtauc_tot = np.maximum(np.asarray(dtauc_tot, dtype=float), 0.0)

        # truncate optical depths
        threshold_od = 1e-8  # threshold at which to truncate optical depths
        significant_layers = np.flatnonzero(dtauc_tot > threshold_od)
        idx = int(significant_layers[0]) if significant_layers.size else 0

        utau = _calculate_output_utau(
            inp.values["out_fmt"],
            output_values,
            altitude_rows,
            dtauc_tot,
            idx,
        )
        model_DISORT.set_utau(utau)

        dtauc_tot = dtauc_tot[idx:]
        tau_g = tau_g[idx:]
        tau_R = tau_R[idx:]
        tau_p = tau_p[idx:]
        w_p = w_p[idx:]
        track_lyr_local = track_lyr[idx:]

        # set layer optical depths
        model_DISORT.set_dtauc_manually(dtauc=dtauc_tot)
        model_DISORT.disort_input["dtauc"][model_DISORT.disort_input["dtauc"] < 0] = 0

        # set maxcly based on current length of tau_g
        model_DISORT.set_maxcly(len(dtauc_tot))

        ## set some maxcly-dependent DISORT input variables

        # set and truncate temper
        model_DISORT.set_temper_from_rfm(model_RFM)
        model_DISORT.disort_input["temper"] = model_DISORT.disort_input["temper"][idx:]

        
        #set bottom boundary temperature
        if "btemp" in inp.values:
            model_DISORT.set_btemp(inp.values["btemp"])
        else:
            model_DISORT.set_btemp(model_DISORT.disort_input["temper"][-1])
        
        # set top boundary temperature    
        if "ttemp" in inp.values:
            model_DISORT.set_ttemp(inp.values["ttemp"])
        else:
            model_DISORT.set_ttemp(model_DISORT.disort_input["temper"][0])
        
        model_DISORT.set_h_lyr(
            np.zeros(shape=(model_DISORT.disort_input["maxcly"] + 1))
        )

        # set single scatter albedo
        model_DISORT.set_ssalb(tau_g=tau_g, tau_R=tau_R, tau_p=tau_p, w_p=w_p)

        particle_moments = {}
        for lyr in particle_lyrs:
            layer_tau = particle_block[lyr]["particle_optical_depth"][block_index]
            layer_ssalb = particle_block[lyr]["single_scattering_albedo"][block_index]
            if layer_tau * layer_ssalb < threshold_od:
                skipped = skipped_scattering[lyr]
                skipped["count"] += 1
                skipped["first"] = wvnm if skipped["first"] is None else skipped["first"]
                skipped["last"] = wvnm
            else:
                coefficients = particle_block[lyr][
                    "normalized_legendre_coefficients"
                ][block_index].copy()
                if (
                    abs(coefficients[0] - 1.0)
                    > layer.NORMALIZED_MOMENT_TOLERANCE
                ):
                    raise RuntimeError(
                        "Something is wrong with the phase function. The first "
                        f"coefficient is {coefficients[0]}, but should be 1.0. "
                        "Try increasing number of quadrature points."
                    )
                coefficients[0] = 1.0
                particle_moments[track_lyr_local.index(lyr)] = coefficients

        model_DISORT.set_mixed_pmom(
            tau_R=tau_R,
            w_p=w_p,
            tau_p=tau_p,
            particle_moments=particle_moments,
            prec=inp.values["disort_precision"],
        )
        #    print(model_DISORT.disort_input["pmom"]).shape

        # set wavenumber range for DISORT (for Planck function)
        model_DISORT.set_wvnm_range(wvnm - 0.5, wvnm + 0.5)

        _set_spectral_boundary_inputs(
            model_DISORT,
            wvl_idx,
            surface_albedo,
            solar_spc,
            sun=inp.values["sun"],
        )
        #    model_DISORT.set_fbeam(0.1)

        # run disort input tests
        #        model_DISORT.test_disort_input_format()
        #        model_DISORT.test_disort_input_integrity()
        # These tests are currently disabled, pending review. It turns out they test
        # inputs, including arrays, element by element, taking about 25% of the total
        # runtime.
        # TODO: options:
        # 1. Remove the tests, potentially unsafe, or rather may be less explanatory
        # when the run fails.
        # 2. Split the tests (i.e. test elements that do not change between runs
        # separately and then test within the loop only the elements that change from
        # iteration to iteration.
        # 3. Move the test behind a debug flag (add the flag to the driver table).
        # 4. Replace element-wise checks with vectorized dtype assertions
        # (e.g. for dtauc do np.asarray(dtauc).dtype / np.issubdtype...

        # call DISOBRDF
        #    model_DISORT.run_disobrdf(prec=inp.values["disort_precision"],
        #                              debug=False,
        #                              brdf_type=2, # Cox-Munk
        #                              brdf_arg=[1,1.34,False,0], # wind speed, water refractive index, do_shadow
        #                              nmug=200) # number of quadrature angles
        # run disort
        disort_result = model_DISORT.run_disort(
            prec=inp.values["disort_precision"],
            adjust_maxcmu=inp.values["adjust_maxcmu"],
        )

        model_SRFM.store_disort_result(disort_result, wvl_idx)

        # print current status
    #    print(model_DISORT.status)
    print("Main DISORT loop finished.")
    for layer_name, skipped in skipped_scattering.items():
        if skipped["count"]:
            warnings.warn(
                f"Particle layer {layer_name!r} was below scattering optical-depth threshold "
                f"{threshold_od:g} at {skipped['count']} wavenumbers from "
                f"{skipped['first']:g} to {skipped['last']:g} cm-1; particle "
                "scattering was omitted there.",
                stacklevel=2,
            )

    # convolve final radiance spectrum with iasi instrument line shape
    if inp.values["convolve_iasi"] == True:
        model_SRFM.convolve_with_iasi(inp.values["iasi_ils"])

    # interpolate resulting bbt and radiances to final grid
    model_SRFM.interp(fin_grid)

    if "bbt" in requested_outputs:
        model_SRFM.calc_bbt()
    if "uu" not in requested_outputs and hasattr(model_SRFM, "uu"):
        delattr(model_SRFM, "uu")
        model_SRFM.retained_outputs = frozenset(
            set(model_SRFM.retained_outputs) - {"uu"}
        )

    ########################################################################################
    # (optional) save spectrum to file(s)
    ########################################################################################
    if inp.values["out_mode"] == "txt":
        if "bbt" in requested_outputs:
            _write_spectral_text(
                os.path.join(inp.values["results_fldr"], "bbt.txt"),
                model_SRFM.wvnm,
                model_SRFM.bbt,
                "Brightness temperature (K)",
            )

        if "uu" in requested_outputs:
            _write_spectral_text(
                os.path.join(inp.values["results_fldr"], "rad.txt"),
                model_SRFM.wvnm,
                model_SRFM.uu,
                "Radiance (W m-2 sr-1 cm)",
            )

    elif inp.values["out_mode"] == "netcdf":
        out_nm = os.path.join(
            inp.values["results_fldr"], inp.values.get("out_fname") or "srfm.nc"
        )
        _write_srfm_netcdf(
            out_nm,
            model_SRFM,
            inp.values["retain_outputs"],
            effective_params,
            RFM_wvnm,
            particle_lyrs,
            grey_body_layers=gbc_lyrs,
            surface_albedo=surface_albedo,
            solar_spectral_irradiance=solar_spc,
            nmom=nmom,
            scattering_block_size=scattering_block_size,
        )

    elif inp.values["out_mode"] == None:
        pass

    if inp.values["base_plots"] == True:
        _plot_retained_spectral_outputs(
            model_SRFM,
            inp.values["results_fldr"],
            inp.values["retain_outputs"],
            inp.values["show_plots"],
        )

    return model_SRFM
