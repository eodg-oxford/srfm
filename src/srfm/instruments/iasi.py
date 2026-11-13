"""IASI-specific processing for :class:`srfm.forward_model.SRFM` results."""

from __future__ import annotations

import os
import warnings
from collections.abc import Mapping
from contextlib import contextmanager
from importlib.resources import as_file, files
from pathlib import Path
from typing import Any, Iterator

import numpy as np
from scipy.signal import convolve

from .. import utilities
from ..forward_model import SRFM


IASI_CHANNEL_COUNT = 8461
IASI_WAVENUMBER_MINIMUM = 645.0
IASI_WAVENUMBER_MAXIMUM = 2760.0
IASI_CHANNEL_SPACING = 0.25
IASI_CONVOLUTION_MARGIN = 2.0
IASI_INSTRUMENT_ALIASES = {
    "m01": "M01",
    "b": "M01",
    "iasi-b": "M01",
    "m02": "M02",
    "a": "M02",
    "iasi-a": "M02",
    "m03": "M03",
    "c": "M03",
    "iasi-c": "M03",
}


def normalize_iasi_instrument(iasi_instrument: str) -> str:
    """Return the data-product identifier for an IASI instrument alias.

    Args:
        iasi_instrument: Canonical identifier or IASI-A, IASI-B, or IASI-C
            alias. Matching ignores surrounding whitespace and letter case.

    Returns:
        The canonical ``"M01"``, ``"M02"``, or ``"M03"`` identifier.

    Raises:
        TypeError: If ``iasi_instrument`` is not a string.
        ValueError: If the string is not a recognized IASI identifier.
    """
    if not isinstance(iasi_instrument, str):
        raise TypeError("IASI instrument must be specified as a string.")
    normalized_alias = iasi_instrument.strip().lower()
    try:
        return IASI_INSTRUMENT_ALIASES[normalized_alias]
    except KeyError as exc:
        raise ValueError(
            "IASI instrument must be M01, M02, M03, a, b, c, "
            "IASI-A, IASI-B, or IASI-C."
        ) from exc


def apply_banded_factor(
    factor_banded: np.ndarray, independent_noise: np.ndarray
) -> np.ndarray:
    """Multiply independent noise by a compact lower-banded factor.

    Args:
        factor_banded: Compact lower-banded factor. Row zero contains the
            diagonal and each subsequent row contains one lower diagonal.
        independent_noise: Independent standard-normal samples, with one
            value for every spectral channel.

    Returns:
        Correlated noise equal to ``L @ independent_noise``.

    Raises:
        ValueError: If the factor is not two-dimensional or the independent
            noise does not contain exactly one value per channel.
    """
    factor_banded = np.asarray(factor_banded)
    independent_noise = np.asarray(independent_noise)

    if factor_banded.ndim != 2:
        raise ValueError("Banded factor must be a two-dimensional array")
    if independent_noise.shape != (factor_banded.shape[1],):
        raise ValueError(
            f"Independent noise must have shape ({factor_banded.shape[1]},)"
        )

    channel_count = factor_banded.shape[1]
    correlated_noise = factor_banded[0] * independent_noise

    for lag in range(1, factor_banded.shape[0]):
        correlated_noise[lag:] += (
            factor_banded[lag, : channel_count - lag]
            * independent_noise[: channel_count - lag]
        )

    return correlated_noise


class IASI:
    """Apply IASI instrument processing to an existing SRFM result.

    Processing methods deliberately mutate the supplied result.  This keeps
    memory use bounded for spectra containing many millions of samples and
    makes the ownership contract explicit to callers.

    Args:
        ils_file: Optional IASI instrument-line-shape file in RFM format. When
            omitted, the copy distributed with SRFM is used.
        iasi_instrument: IASI instrument identifier or accepted alias.
        iasi_ifov: Instrument field of view, from one to four. The default is
            one.
    """

    name = "iasi"

    def __init__(
        self,
        ils_file: str | os.PathLike[str] | None = None,
        iasi_instrument: str = "M03",
        iasi_ifov: int = 1,
    ):
        self.ils_file = None if ils_file is None else Path(ils_file)
        self.iasi_instrument = normalize_iasi_instrument(iasi_instrument)
        self.iasi_ifov = iasi_ifov

    @contextmanager
    def _resolved_ils_path(self) -> Iterator[Path]:
        """Yield an explicit ILS path or the packaged default resource."""
        if self.ils_file is not None:
            yield self.ils_file
            return

        resource = files("srfm").joinpath("RFM", "rfm_files", "iasi.ils")
        if not resource.is_file():
            raise FileNotFoundError(
                "The packaged IASI instrument line shape is unavailable."
            )
        with as_file(resource) as packaged_path:
            yield packaged_path

    @contextmanager
    def _resolved_noise_path(self) -> Iterator[Path]:
        """Yield the installed IASI noise-model resource as a file path."""
        resource = files("srfm").joinpath("data", "iasi_noise.npz")
        if not resource.is_file():
            raise FileNotFoundError(
                "The packaged IASI noise model 'iasi_noise.npz' is unavailable."
            )
        with as_file(resource) as packaged_path:
            yield packaged_path

    @property
    def ils_source(self) -> str:
        """Return a stable label for metadata and diagnostics."""
        if self.ils_file is None:
            return "package:srfm/RFM/rfm_files/iasi.ils"
        return os.fspath(self.ils_file)

    @classmethod
    def from_inputs(cls, values: Mapping[str, Any]) -> IASI:
        """Construct IASI from a validated runner input mapping.

        The generic and OXHARP runners use ``iasi_ils`` while the tailored
        IASI runner retains its established ``ils`` input name.  Keeping that
        translation here prevents the runners and generic selector from
        depending on IASI-specific parameter names.
        """
        ils_file = values.get("iasi_ils") or values.get("ils")
        if values.get("add_noise") == "iasi":
            iasi_instrument = values.get("iasi_instrument", "M03")
            iasi_ifov = values.get("iasi_ifov", 1)
        else:
            iasi_instrument = "M03"
            iasi_ifov = 1
        return cls(ils_file, iasi_instrument, iasi_ifov)

    def output_grid(self, values: Mapping[str, Any]) -> np.ndarray:
        """Derive the IASI output channels from the computational range.

        Two inverse centimetres are removed from both computational boundaries
        so the returned channels do not use the convolution's affected edge
        region. Configured final-grid values are ignored because IASI output is
        always sampled at the instrument's channel locations.

        Args:
            values: Validated runner inputs containing the computational
                spectral boundaries and units.

        Returns:
            One-dimensional array of IASI channel wavenumbers in inverse
            centimetres.

        Raises:
            ValueError: If the computational range is not expressed in inverse
                centimetres, extends beyond the supported IASI range including
                its margins, is narrower than four inverse centimetres, or
                contains no IASI output channel after removing the margins.
        """
        if values.get("spc_units") != "cm-1":
            raise ValueError(
                "IASI convolution requires spc_units='cm-1' because its "
                "two-sided margin and instrument line shape are specified in "
                "inverse centimetres."
            )

        computational_lower = float(values["spc_wvnmlo"])
        computational_upper = float(values["spc_wvnmhi"])
        computational_width = computational_upper - computational_lower
        if computational_width < 2.0 * IASI_CONVOLUTION_MARGIN:
            raise ValueError(
                "IASI convolution requires spc_wvnmhi - spc_wvnmlo to be at "
                "least 4 cm-1 so a 2 cm-1 margin can be removed from each side."
            )

        supported_computational_lower = (
            IASI_WAVENUMBER_MINIMUM - IASI_CONVOLUTION_MARGIN
        )
        supported_computational_upper = (
            IASI_WAVENUMBER_MAXIMUM + IASI_CONVOLUTION_MARGIN
        )
        if (
            computational_lower < supported_computational_lower
            or computational_upper > supported_computational_upper
        ):
            raise ValueError(
                "IASI convolution requires the computational spectral range to "
                "lie within 643 to 2762 cm-1, including the 2 cm-1 convolution "
                "margin."
            )

        requested_output_lower = (
            computational_lower + IASI_CONVOLUTION_MARGIN
        )
        requested_output_upper = (
            computational_upper - IASI_CONVOLUTION_MARGIN
        )
        complete_iasi_grid = IASI_WAVENUMBER_MINIMUM + np.arange(
            IASI_CHANNEL_COUNT, dtype=float
        ) * IASI_CHANNEL_SPACING
        boundary_tolerance = 1e-10
        selected_channels = complete_iasi_grid[
            (complete_iasi_grid >= requested_output_lower - boundary_tolerance)
            & (complete_iasi_grid <= requested_output_upper + boundary_tolerance)
        ]
        if selected_channels.size == 0:
            raise ValueError(
                "The computational spectral range contains no IASI channel after "
                "removing the 2 cm-1 convolution margin from each side."
            )

        ignored_final_fields = tuple(
            field_name
            for field_name in ("fin_wvnmlo", "fin_wvnmhi", "fin_res")
            if values.get(field_name) is not None
        )
        if ignored_final_fields:
            warnings.warn(
                "convolve='iasi' ignores "
                + ", ".join(ignored_final_fields)
                + "; the actual output contains IASI channels from "
                f"{selected_channels[0]:g} to {selected_channels[-1]:g} cm-1, "
                "derived from spc_wvnmlo and spc_wvnmhi after removing a "
                "2 cm-1 margin from each side.",
                UserWarning,
                stacklevel=2,
            )
        return selected_channels

    def convolve(self, result: SRFM) -> SRFM:
        """Convolve every retained radiance spectrum in place.

        One spectrum is convolved at a time and written back into the existing
        ``uu`` array.  The method therefore avoids allocating a second array the
        size of the complete multidimensional radiance result.

        Args:
            result: SRFM result with a regular wavenumber grid and ``uu`` whose
                first dimension is spectral.

        Returns:
            The same ``result`` instance supplied by the caller.

        Raises:
            TypeError: If ``result`` is not an SRFM instance.
            ValueError: If its spectral grid, radiance shape, or ILS is invalid.
        """
        if not isinstance(result, SRFM):
            raise TypeError("IASI.convolve expects an SRFM result.")
        if not hasattr(result, "uu"):
            raise ValueError("IASI convolution requires retained radiance 'uu'.")
        if not hasattr(result, "wvnm"):
            raise ValueError("IASI convolution requires a wavenumber grid.")

        wavenumber = np.asarray(result.wvnm, dtype=float)
        radiance = np.asarray(result.uu)
        if wavenumber.ndim != 1 or wavenumber.size < 2:
            raise ValueError("The wavenumber grid must contain at least two points.")
        if not np.all(np.isfinite(wavenumber)):
            raise ValueError("The wavenumber grid must contain only finite values.")
        if radiance.ndim < 1 or radiance.shape[0] != wavenumber.size:
            raise ValueError(
                "The first radiance dimension must match the wavenumber grid."
            )
        if not np.issubdtype(radiance.dtype, np.floating):
            raise TypeError("IASI convolution requires floating-point radiance data.")

        spacing = np.diff(wavenumber)
        if not np.allclose(spacing, spacing[0], rtol=1e-10, atol=1e-12):
            raise ValueError("Wavenumber grid is not regular.")
        # Retain the historical eight-decimal ILS resampling convention so the
        # architectural extraction does not alter existing IASI calculations.
        resolution = float(np.round(spacing[0], decimals=8))
        if resolution <= 0:
            raise ValueError(
                "Wavenumber grid must be strictly increasing with resolution "
                "of at least 1e-8 cm-1."
            )

        with self._resolved_ils_path() as ils_path:
            ils_x, ils_y, ils_low, ils_high = utilities.read_ils(
                os.fspath(ils_path)
            )
        if not (
            np.all(np.isfinite(ils_x))
            and np.all(np.isfinite(ils_y))
            and len(ils_x) == len(ils_y)
        ):
            raise ValueError("Instrument line shape must contain finite paired values.")

        point_count = int(np.floor((ils_high - ils_low) / resolution)) + 1
        if point_count < 1:
            raise ValueError("Instrument line shape does not cover a usable interval.")
        response_grid = ils_low + np.arange(point_count) * resolution
        response = np.interp(response_grid, ils_x, ils_y)
        normalization = np.sum(response)
        if not np.isfinite(normalization) or np.isclose(normalization, 0.0):
            raise ValueError(
                "Instrument line shape has zero or non-finite normalization."
            )

        # Check before reshaping so an unusual layout cannot cause NumPy to
        # allocate a complete copy merely to discover that it is unsuitable.
        if not radiance.flags.c_contiguous:
            raise ValueError(
                "Radiance layout cannot be processed in place without a full copy."
            )
        flattened = radiance.reshape(wavenumber.size, -1)
        for spectrum_index in range(flattened.shape[1]):
            convolved = (
                convolve(flattened[:, spectrum_index], response, mode="same")
                / normalization
            )
            flattened[:, spectrum_index] = convolved

        # Discard an existing derived array before allocating its replacement.
        if hasattr(result, "bbt"):
            delattr(result, "bbt")
            result.calc_bbt()

        result.instrument_name = self.name
        result.spectral_representation = "iasi_ils_convolved_spectral_radiance"
        result.spectral_coordinate_role = "wavenumber"
        result.noise_applied = False
        history = tuple(getattr(result, "processing_history", ()))
        result.processing_history = history + (
            {
                "operation": "convolution",
                "instrument": self.name,
                "ils_file": self.ils_source,
            },
        )
        return result

    def add_noise(
        self,
        result: SRFM,
        rng: np.random.Generator | None = None,
    ) -> SRFM:
        """Add correlated IASI radiance noise to every spectrum in place.

        The packaged lower-banded factor is applied separately to every
        combination of polar angle, output level, and azimuthal angle. Only
        two channel-length temporary arrays are needed for each spectrum.

        Args:
            result: Convolved SRFM result on a contiguous part of the IASI
                channel grid.
            rng: Optional NumPy random generator. A fresh default generator is
                created when omitted.

        Returns:
            The same ``result`` instance supplied by the caller.

        Raises:
            TypeError: If ``result`` is not an SRFM instance, its radiance is
                not floating point, or ``rng`` is not a NumPy generator.
            ValueError: If radiance, spectral grid, instrument, or IFOV
                configuration is incompatible with the packaged noise model.
            FileNotFoundError: If the packaged noise model is unavailable.
        """
        if not isinstance(result, SRFM):
            raise TypeError("IASI.add_noise expects an SRFM result.")
        if not hasattr(result, "uu"):
            raise ValueError(
                "IASI noise application requires retained radiance 'uu'."
            )
        if not hasattr(result, "wvnm"):
            raise ValueError("IASI noise application requires a wavenumber grid.")
        if getattr(result, "spectral_representation", None) != (
            "iasi_ils_convolved_spectral_radiance"
        ):
            raise ValueError(
                "IASI noise application requires an IASI-convolved spectrum; "
                "set convolve='iasi' when add_noise='iasi'."
            )
        if (
            type(self.iasi_ifov) is not int
            or self.iasi_ifov not in {1, 2, 3, 4}
        ):
            raise ValueError("IASI IFOV must be one of 1, 2, 3, or 4 to add noise.")
        if rng is None:
            rng = np.random.default_rng()
        elif not isinstance(rng, np.random.Generator):
            raise TypeError("rng must be a numpy.random.Generator instance.")

        radiance = np.asarray(result.uu)
        result_wavenumber_cm_inverse = np.asarray(result.wvnm, dtype=float)
        if radiance.ndim < 1 or radiance.shape[0] != result_wavenumber_cm_inverse.size:
            raise ValueError(
                "The first radiance dimension must match the IASI output grid."
            )
        if not np.issubdtype(radiance.dtype, np.floating):
            raise TypeError(
                "IASI noise application requires floating-point radiance."
            )
        if (
            result_wavenumber_cm_inverse.ndim != 1
            or result_wavenumber_cm_inverse.size == 0
        ):
            raise ValueError(
                "IASI noise application requires a non-empty one-dimensional "
                "wavenumber grid."
            )

        noise_key = f"{self.iasi_instrument}_{self.iasi_ifov}"
        with self._resolved_noise_path() as noise_path:
            with np.load(noise_path) as noise_data:
                if noise_key not in noise_data:
                    raise ValueError(
                        f"The packaged IASI noise model has no factor for {noise_key}."
                    )
                factor_banded = noise_data[noise_key]
                noise_wavenumber_cm_inverse = (
                    noise_data["wavenumber_m_inverse"] / 100.0
                )

        if factor_banded.shape[1] != IASI_CHANNEL_COUNT:
            raise ValueError(
                f"IASI noise factor {noise_key} has {factor_banded.shape[1]} "
                "channels; expected 8461."
            )
        channel_indices = np.searchsorted(
            noise_wavenumber_cm_inverse, result_wavenumber_cm_inverse
        )
        indices_are_in_range = np.all(channel_indices < IASI_CHANNEL_COUNT)
        if indices_are_in_range:
            matched_wavenumber = noise_wavenumber_cm_inverse[channel_indices]
            channels_match = np.allclose(
                result_wavenumber_cm_inverse,
                matched_wavenumber,
                rtol=0.0,
                atol=1e-10,
            )
        else:
            channels_match = False
        channels_are_contiguous = (
            channel_indices.size == 1
            or np.all(np.diff(channel_indices) == 1)
        )
        if not channels_match or not channels_are_contiguous:
            raise ValueError(
                "The SRFM wavenumber grid must be a contiguous subset of the "
                "packaged IASI grid from 645 to 2760 cm-1 at 0.25 cm-1 spacing."
            )

        spectrum_dimensions = radiance.shape[1:]
        if spectrum_dimensions:
            spectrum_indices = np.ndindex(spectrum_dimensions)
        else:
            spectrum_indices = [()]
        for spectrum_index in spectrum_indices:
            independent_noise = rng.standard_normal(IASI_CHANNEL_COUNT)
            correlated_noise = apply_banded_factor(
                factor_banded, independent_noise
            )
            radiance[(slice(None),) + spectrum_index] += correlated_noise[
                channel_indices
            ]

        if hasattr(result, "bbt"):
            delattr(result, "bbt")
            result.calc_bbt()

        result.instrument_name = self.name
        result.iasi_instrument = self.iasi_instrument
        result.iasi_ifov = self.iasi_ifov
        if result.spectral_representation == "iasi_ils_convolved_spectral_radiance":
            result.spectral_representation = (
                "iasi_ils_convolved_noisy_spectral_radiance"
            )
        else:
            result.spectral_representation = "iasi_noisy_spectral_radiance"
        result.noise_applied = True
        history = tuple(getattr(result, "processing_history", ()))
        result.processing_history = history + (
            {
                "operation": "noise",
                "instrument": self.name,
                "iasi_instrument": self.iasi_instrument,
                "iasi_ifov": self.iasi_ifov,
                "noise_factor": noise_key,
            },
        )
        return result
