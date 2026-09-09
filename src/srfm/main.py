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


def _interpolate_scattering_block(scattering_layers, wavelengths):
    """Interpolate every particle layer for one bounded spectral block.

    Args:
        scattering_layers (Mapping[str, MieLayer]): Coarse-grid particle layers.
        wavelengths (array-like): Wavelengths for the current block.

    Returns:
        dict[str, dict[str, numpy.ndarray]]: Extinction, optical depth, albedo,
        and Legendre coefficients for each named layer.
    """
    block = {}
    for layer_name, scattering_layer in scattering_layers.items():
        properties = scattering_layer.interpolate_optical_properties(wavelengths)
        properties["tau"] = (
            properties["beta_ext"]
            * 1e3
            * (scattering_layer.alt_upp - scattering_layer.alt_low)
        )
        block[layer_name] = properties
    return block


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


def _write_srfm_netcdf(
    filename,
    model,
    retained_outputs,
    effective_params,
    wavelengths,
    scattering_layers,
):
    """Write one complete SRFM NetCDF file.

    The file contains every explicitly retained SRFM result, shared output
    coordinates, and the existing scattering-layer optical-depth metadata.

    Args:
        filename (path-like): Destination NetCDF filename.
        model (SRFM): Completed model containing retained result arrays.
        retained_outputs (collection[str]): User-selected result names.
        effective_params (Mapping): Effective validated run configuration.
        wavelengths (array-like): Computational wavelength grid.
        scattering_layers (Mapping[str, MieLayer]): Configured scattering layers.

    Returns:
        None: The NetCDF file is written and closed.
    """
    wavelengths = np.asarray(wavelengths)
    layer_names = list(scattering_layers)
    serialized_params = copy.deepcopy(effective_params)
    serialized_params["driver_inputs"]["spectral"] = str(
        serialized_params["driver_inputs"]["spectral"]
    )

    with Dataset(filename, "w", format="NETCDF4") as nc_file:
        nc_file.description = "SRFM output."
        nc_file.history = (
            f"Created {datetime.datetime.now().strftime('%Y-%m-%d')}"
        )
        nc_file.srfm_params = json.dumps(
            serialized_params,
            default=utilities.json_handler,
        )

        nc_file.createDimension("wavenumber_op", wavelengths.size)
        nc_file.createDimension("layer", len(layer_names))
        _create_netcdf_spectral_dimensions(nc_file, model)
        _write_retained_netcdf_outputs(nc_file, model, retained_outputs)

        layer_variable = nc_file.createVariable("layer_names", str, ("layer",))
        layer_variable[:] = np.asarray(layer_names, dtype=object)

        optical_depth = np.zeros((len(layer_names), wavelengths.size))
        for index, layer_name in enumerate(layer_names):
            scattering_layer = scattering_layers[layer_name]
            beta_ext = scattering_layer.interpolate_optical_properties(
                wavelengths, include_legendre=False
            )["beta_ext"]
            optical_depth[index, :] = (
                beta_ext
                * 1e3
                * (scattering_layer.alt_upp - scattering_layer.alt_low)
            )

        tau = nc_file.createVariable(
            "tau",
            "f8",
            ("layer", "wavenumber_op"),
            zlib=True,
            complevel=4,
        )
        tau.long_name = "Optical depth per layer and wavenumber"
        tau[:] = optical_depth


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
    # Assign some variables:
    ########################################################################################
    # Keep all run-generated files outside the installed package tree.
    os.makedirs(inp.values["results_fldr"], exist_ok=True)

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

    rfm_grid_fname = rfm_functions.construct_rfm_grid_file(
        RFM_wvnm, filename="grid.spc", rfm_fldr=inp.values["results_fldr"]
    )
    ########################################################################################
    # define an atmospheric scattering layers
    ########################################################################################
    scat_lyrs = (
        {}
    )  # dictionary of scattering layers, key - layer name, value - Layer() object

    # define Layer properties
    scat_lyrs_inputs = {}

    if "scat_lyrs_inputs" in inp.values.keys():

        # calculate MieLayer optical properties
        for lyr in inp.values["scat_lyrs_inputs"].keys():
            scat_lyrs_inputs[lyr] = (
                inp.values["scat_lyrs_inputs"][lyr]
                | inp.values["scat_lyrs_inputs"][lyr]
            )
            # note: the pipe "|" here creates a shallow copy, i.e.
            # inp.values["scat_lyers_inputs"][lyr] and scat_layers_inputs][lyr]
            # are now different objects in memory
            scat_lyrs[lyr] = layer.MieLayer()
            scat_lyrs[lyr].set_input_from_dict(
                scat_lyrs_inputs[lyr]
            )  # sets input for scattering layer
            scat_lyrs[
                lyr
            ].calculate_op()  # calculates layer optical properties, may run in parallel

    # add output from optical properties calculation
    for lyr in scat_lyrs.keys():
        scat_lyrs[lyr].add_op_calc_output()
        if not inp.values.get("retain_phase_functions", False):
            scat_lyrs[lyr].discard_phase_function()

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

    effective_params = copy.deepcopy(inp.values)
    effective_params["scat_lyrs_inputs"] = {
      lyr: {
          attr: getattr(scat_lyrs[lyr], attr)
          for attr in layer_attrs
          if hasattr(scat_lyrs[lyr], attr)
      }
      for lyr in scat_lyrs
    }

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
    for lyr in scat_lyrs.keys():
        if (scat_lyrs[lyr].legendre_coefficient.shape[1] - 1) > nmom:
            nmom = scat_lyrs[lyr].legendre_coefficient.shape[1] - 1

    model_DISORT.set_maxcmu(inp.values["maxcmu"])

    model_DISORT.set_maxmom(nmom)
    if nmom < model_DISORT.disort_input["maxcmu"]:
        model_DISORT.set_maxmom(model_DISORT.disort_input["maxcmu"])

    model_DISORT.set_maxumu(inp.values["maxumu"])
    model_DISORT.set_maxphi(inp.values["maxphi"])
    # maxulv must match the resolved output geometry. The driver-table value is
    # retained for compatibility with other runners but is derived here.
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
    model_DISORT.set_albedo(inp.values["albedo"])

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
    if inp.values["sun"] == True:
        # load solar spectrum from file
        solar_spc, solar_spc_wvnm = utilities.load_solar_spectrum_Gueymard20018()
        solar_spc = np.interp(
            RFM_wvnm, solar_spc_wvnm[::-1], solar_spc[::-1]
        )  # interpolate to RFM_wvnm (the calculation grid)

        # scale with year day (different Sun-Earth distance throughout the year
        # the original spectrum is for 1 AU
        solar_spc = utilities.scale_solar_spectrum(solar_spc, year_day)

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
    scattering_block = {}
    scattering_block_start = 0
    skipped_scattering = {
        layer_name: {"count": 0, "first": None, "last": None}
        for layer_name in scat_lyrs
    }
    for wvl_idx, (wvnm, wvl, tau_g) in enumerate(
        zip(RFM_wvnm, wvls, model_RFM.rfm_output.differential_tau)
    ):
        if scat_lyrs and wvl_idx % scattering_block_size == 0:
            scattering_block_start = wvl_idx
            block_stop = min(wvl_idx + scattering_block_size, len(wvls))
            scattering_block = _interpolate_scattering_block(
                scat_lyrs, wvls[wvl_idx:block_stop]
            )
        block_index = wvl_idx - scattering_block_start

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

        tau_g = np.asarray(tau_g, dtype=float)

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
        for lyr in scat_lyrs.keys():
            tau_p[track_lyr.index(lyr)] = scattering_block[lyr]["tau"][block_index]

        # particle layer single scatter albedo
        w_p = np.zeros(shape=(len(tau_g)))
        for lyr in scat_lyrs.keys():
            w_p[track_lyr.index(lyr)] = scattering_block[lyr]["ssalb"][block_index]

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
        for lyr in scat_lyrs.keys():
            layer_tau = scattering_block[lyr]["tau"][block_index]
            if layer_tau < threshold_od:
                skipped = skipped_scattering[lyr]
                skipped["count"] += 1
                skipped["first"] = wvnm if skipped["first"] is None else skipped["first"]
                skipped["last"] = wvnm
            else:
                coefficients = scattering_block[lyr]["legendre_coefficient"][
                    block_index
                ].copy()
                Legendre_precision = 1 / coefficients[0]
                if abs(Legendre_precision - 1) > 1e-5:
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

        # set incoming beam of (solar) radiation
        if inp.values["sun"] == True:
            model_DISORT.set_fbeam(solar_spc[wvl_idx])
        else:
            model_DISORT.set_fbeam(0)
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
                f"Scattering layer {layer_name!r} was below optical-depth threshold "
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
            wvls,
            scat_lyrs,
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
