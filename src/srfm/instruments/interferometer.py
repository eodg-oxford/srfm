"""Channel-domain processing shared by CrIS and the two HIRAS generations.

The spectral integration is performed on radiances.  In particular, no Fourier
transform of the calculated spectrum and no dense channel-by-input matrix are
needed.  The response and noise data bundled with this module are nominal
product models; their provenance is recorded in the result.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import cumulative_trapezoid

from ..forward_model import SRFM


CHANNEL_STEP = 0.625
HAMMING = np.array((0.23, 0.54, 0.23))


@dataclass(frozen=True)
class Band:
    name: str
    first: float
    count: int
    response: str = "sinc"

    @property
    def channels(self) -> np.ndarray:
        return self.first + CHANNEL_STEP * np.arange(self.count)


def _key(value: float) -> int:
    """Use exact 1/1000 cm-1 integer keys for overlapping band coordinates."""
    return round(float(value) * 1000)


@lru_cache(maxsize=12)
def _response_for(band: Band, platform: str) -> tuple[np.ndarray | None, np.ndarray | None, float]:
    if band.response == "sinc":
        return None, None, 30.0
    resource = files("srfm").joinpath("data", "hiras2_srf.npz")
    with np.load(resource) as data:
        prefix = "FY3E_LW" if platform == "FY3E" else platform
        offsets = np.asarray(data[prefix + "_offset"], dtype=float)
        response = np.asarray(data[prefix + "_response"], dtype=float)
    if offsets.ndim != 1 or not np.all(np.diff(offsets) > 0):
        raise ValueError(f"Invalid {platform} {band.name} spectral response offsets.")
    if not np.all(np.isfinite(response)):
        raise ValueError(f"Invalid {platform} {band.name} spectral response values.")
    return offsets, response, float(max(abs(offsets[0]), abs(offsets[-1])))


def _response_values(
    offsets: np.ndarray,
    table_offsets: np.ndarray | None,
    table_response: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the signed response and its slope for moment-aware quadrature."""
    if table_offsets is None:
        value = np.sinc(offsets / CHANNEL_STEP)
        delta = 1e-4
        slope = (
            np.sinc((offsets + delta) / CHANNEL_STEP)
            - np.sinc((offsets - delta) / CHANNEL_STEP)
        ) / (2 * delta)
    else:
        value = np.interp(offsets, table_offsets, table_response)
        table_slope = np.gradient(table_response, table_offsets)
        slope = np.interp(offsets, table_offsets, table_slope)
    return value, slope


def _integrate_channels(
    x: np.ndarray,
    y: np.ndarray,
    centers: np.ndarray,
    table_offsets: np.ndarray | None,
    table_response: np.ndarray | None,
    support: float,
) -> np.ndarray:
    """Integrate one band for all retained spatial spectra using bounded memory."""
    left = float(centers[0] - support)
    right = float(centers[-1] + support)
    if x[0] > left + 1e-9 or x[-1] < right - 1e-9:
        raise ValueError(
            f"Convolution needs high-resolution coverage {left:g} to "
            f"{right:g} cm-1; supplied range is {x[0]:g} to {x[-1]:g} cm-1."
        )
    low = max(0, int(np.searchsorted(x, left, side="right")) - 1)
    high = min(x.size, int(np.searchsorted(x, right, side="left")) + 1)
    xx = x[low:high]
    yy = y[low:high].reshape(high - low, -1)
    source_step = float(np.median(np.diff(xx)))
    bin_step = max(source_step, 0.01)
    nbin = int(np.ceil((right - left) / bin_step))
    edges = np.linspace(left, right, nbin + 1)
    mids = (edges[:-1] + edges[1:]) / 2
    widths = np.diff(edges)

    # Store only response weights, not a channel by full-spectrum matrix.
    weights: list[tuple[np.ndarray, np.ndarray, np.ndarray, float]] = []
    for center in centers:
        start = int(np.searchsorted(mids, center - support, side="left"))
        stop = int(np.searchsorted(mids, center + support, side="right"))
        relative = mids[start:stop] - center
        response, slope = _response_values(relative, table_offsets, table_response)
        area = float(np.dot(response, widths[start:stop]))
        if not np.isfinite(area) or abs(area) < 1e-10:
            raise ValueError(f"Response at {center:g} cm-1 has invalid area {area:g}.")
        weights.append((start, stop, response / area, slope / area))

    result = np.empty((centers.size, yy.shape[1]), dtype=float)
    for spectrum in range(yy.shape[1]):
        radiance = yy[:, spectrum]
        cumulative = cumulative_trapezoid(radiance, xx, initial=0)
        cumulative_first = cumulative_trapezoid(radiance * xx, xx, initial=0)
        integral = np.diff(np.interp(edges, xx, cumulative))
        first = np.diff(np.interp(edges, xx, cumulative_first)) - mids * integral
        for row, (start, stop, response, slope) in enumerate(weights):
            result[row, spectrum] = (
                np.dot(response, integral[start:stop])
                + np.dot(slope, first[start:stop])
            )
    return result.reshape((centers.size,) + y.shape[1:])


class Interferometer:
    """Common native-channel convolution and FOV-specific noise handling."""

    name = ""
    default_noise_resource = ""
    fov_count = 0

    def __init__(
        self,
        *,
        platform: str,
        fov: int,
        bands: tuple[Band, ...],
        values: dict[str, Any] | None = None,
        apodization: str | None = None,
        noise_file: str | None = None,
        noise_units: str | None = None,
        noise_seed: int | None = None,
    ) -> None:
        self.platform = platform
        self.fov = fov
        self.bands = bands
        self.values = values
        self.apodization = apodization
        self.noise_file = noise_file
        self.noise_units = noise_units
        self.noise_seed = noise_seed
        self._selected: list[tuple[int, int]] | None = None
        if type(fov) is not int or not 1 <= fov <= self.fov_count:
            raise ValueError(
                f"{self.name} {platform}: FOV must be 1..{self.fov_count}; got {fov!r}."
            )
        if apodization not in (None, "Hamming"):
            raise ValueError(
                f"{self.name} {platform}: apodization must be None or 'Hamming'."
            )

    def _all_unique(self) -> list[tuple[float, int, int]]:
        # Earlier bands own the eight FY-3E/F overlap coordinates.
        owner: dict[int, tuple[float, int, int]] = {}
        for band_number, band in enumerate(self.bands):
            for index, center in enumerate(band.channels):
                if self.apodization == "Hamming" and (
                    index < 2 or index >= band.count - 2
                ):
                    continue
                owner.setdefault(_key(center), (float(center), band_number, index))
        return sorted(owner.values())

    def output_grid(self, values: dict[str, Any]) -> np.ndarray:
        if values.get("spc_units") != "cm-1":
            raise ValueError(
                f"{self.name} {self.platform}: spc_units must be 'cm-1' for "
                "instrument channels and response support."
            )
        resolution = values.get("fin_res")
        if resolution is not None and not np.isclose(resolution, CHANNEL_STEP, rtol=0, atol=1e-10):
            raise ValueError(
                f"{self.name} {self.platform} FSR: fin_res={resolution:g} cm-1 "
                f"is incompatible with its {CHANNEL_STEP:g} cm-1 channel grid; "
                "set fin_res=None and use fin_wvnmlo/hi to select channels."
            )
        if float(values["spc_res"]) > 0.05:
            raise ValueError(
                f"{self.name} {self.platform}: spc_res={values['spc_res']:g} "
                "cm-1 is too coarse for the nominal response; use <=0.05 cm-1 "
                "and check convergence for line-rich spectra."
            )
        lower = values.get("fin_wvnmlo")
        upper = values.get("fin_wvnmhi")
        all_channels = self._all_unique()
        if lower is not None and lower < all_channels[0][0] - 1e-9:
            raise ValueError(
                f"{self.name} {self.platform}: fin_wvnmlo={lower:g} is below "
                f"the first instrument channel {all_channels[0][0]:g} cm-1."
            )
        if upper is not None and upper > all_channels[-1][0] + 1e-9:
            raise ValueError(
                f"{self.name} {self.platform}: fin_wvnmhi={upper:g} is above "
                f"the last instrument channel {all_channels[-1][0]:g} cm-1."
            )
        selected = []
        missing = []
        spc_low = float(values["spc_wvnmlo"])
        spc_high = float(values["spc_wvnmhi"])
        for center, band_number, index in all_channels:
            if lower is not None and center < lower - 1e-9:
                continue
            if upper is not None and center > upper + 1e-9:
                continue
            band = self.bands[band_number]
            if self.apodization == "Hamming" and (index < 2 or index >= band.count - 2):
                continue
            _, _, support = _response_for(band, self.platform)
            margin = CHANNEL_STEP if self.apodization == "Hamming" else 0.0
            if center - support - margin < spc_low - 1e-9 or center + support + margin > spc_high + 1e-9:
                missing.append((center, support + margin))
                continue
            selected.append((center, band_number, index))
        if missing and (lower is not None or upper is not None):
            center, margin = missing[0]
            raise ValueError(
                f"{self.name} {self.platform}: channel {center:g} cm-1 needs "
                f"computational coverage {center-margin:g} to {center+margin:g} "
                f"cm-1 for its response{' and Hamming neighbour' if self.apodization else ''}; "
                f"spc range is {spc_low:g} to {spc_high:g} cm-1."
            )
        if not selected:
            raise ValueError(
                f"{self.name} {self.platform}: no native channels are available "
                f"within the requested range and response margins; spc range "
                f"is {spc_low:g} to {spc_high:g} cm-1."
            )
        self._selected = [(b, i) for _, b, i in selected]
        return np.asarray([c for c, _, _ in selected])

    def convolve(self, result: SRFM) -> SRFM:
        if self.values is None:
            raise ValueError("Standalone convolution requires instrument grid inputs.")
        return self.process_to_grid(result, self.output_grid(self.values))

    def process_to_grid(self, result: SRFM, grid: np.ndarray) -> SRFM:
        """Compute instrument channels, interpolate other fields, and replace radiance."""
        if not isinstance(result, SRFM) or not hasattr(result, "uu"):
            raise ValueError(f"{self.name} convolution needs an SRFM radiance result.")
        if self.values is None or not np.array_equal(grid, self.output_grid(self.values)):
            raise ValueError("Instrument output grid does not match configured channels.")
        source_grid = np.asarray(result.wvnm, dtype=float)
        source = np.asarray(result.uu)
        if (source_grid.ndim != 1 or source.ndim < 1 or
                source.shape[0] != source_grid.size or
                not np.all(np.isfinite(source_grid)) or
                not np.all(np.diff(source_grid) > 0)):
            raise ValueError("Source radiance must have a strictly increasing matching wavenumber grid.")
        computed = np.empty((len(grid),) + source.shape[1:], dtype=float)
        for band_number, band in enumerate(self.bands):
            output_rows = [row for row, (b, _) in enumerate(self._selected) if b == band_number]
            if not output_rows:
                continue
            selected_indices = [self._selected[row][1] for row in output_rows]
            needed = set(selected_indices)
            if self.apodization == "Hamming":
                for index in selected_indices:
                    needed.update((index - 1, index + 1))
            native_indices = np.asarray(sorted(needed), dtype=int)
            tab_x, tab_y, support = _response_for(band, self.platform)
            native = _integrate_channels(
                source_grid, source, band.channels[native_indices], tab_x, tab_y, support
            )
            native_map = {int(index): row for row, index in enumerate(native_indices)}
            for output_row, index in zip(output_rows, selected_indices):
                if self.apodization == "Hamming":
                    computed[output_row] = (
                        HAMMING[0] * native[native_map[index - 1]]
                        + HAMMING[1] * native[native_map[index]]
                        + HAMMING[2] * native[native_map[index + 1]]
                    )
                else:
                    computed[output_row] = native[native_map[index]]
        result.interp(grid)
        result.uu = computed
        if hasattr(result, "bbt"):
            result.calc_bbt()
        result.instrument_name = self.name
        result.instrument_platform = self.platform
        result.instrument_fov = self.fov
        result.instrument_apodization = self.apodization
        result.spectral_representation = f"{self.name}_{'hamming' if self.apodization else 'native'}_channel_radiance"
        result.noise_applied = False
        result.processing_history = tuple(getattr(result, "processing_history", ())) + ({
            "operation": "convolution", "instrument": self.name,
            "platform": self.platform, "fov": self.fov,
            "apodization": self.apodization or "none",
            "quadrature_bin_target_cm-1": 0.01,
            "responses": {
                band.name: (
                    "package:hiras2_srf.npz" if band.response == "table"
                    else "nominal_sinc_MPD_0.8_cm"
                ) for band in self.bands
            },
        },)
        return result

    def _noise_arrays(self) -> tuple[dict[str, np.ndarray], str]:
        if self.noise_file is None:
            if not self.default_noise_resource:
                raise ValueError(
                    f"{self.name} {self.platform}: add_noise requires a matching "
                    "noise_file; no platform-specific noise profile is packaged."
                )
            path = files("srfm").joinpath("data", self.default_noise_resource)
            source = "package:" + self.default_noise_resource
        else:
            path = Path(self.noise_file)
            source = str(path)
        if str(path).lower().endswith(".npz"):
            with np.load(path) as profile:
                data = {key: profile[key] for key in profile.files}
            units = str(data.get("units", ""))
            source_platform = str(data.get("platform", ""))
            if "source" in data:
                source += " (" + str(data["source"]) + ")"
        else:
            from netCDF4 import Dataset
            filename = Path(path).name.upper().replace("-", "")
            source_platform = (
                "J2" if ".J2.CRIS." in filename or ".NOAA21.CRIS." in filename else
                next((p for p in ("FY3D", "FY3E", "FY3F", "FY3H")
                      if filename.startswith(p)), "")
            )
            with Dataset(path) as product:
                product.set_auto_mask(False)
                group = product.groups.get("Data", product)
                group.set_auto_mask(False)
                data = {}
                observed_units = set()
                for band in self.bands:
                    if self.name == "cris":
                        x = group.variables["wnum_" + band.name.lower()][:]
                        variable = group.variables["nedn_" + band.name.lower()]
                        sigma = variable[:]
                    else:
                        x = group.variables["WL_" + band.name][:]
                        variable = group.variables["DS_NEdN_" + band.name]
                        raw = np.asarray(variable[:], dtype=float)
                        fill = getattr(variable, "FillValue", None)
                        valid = np.isfinite(raw) & (raw > 0)
                        if fill is not None:
                            valid &= raw != float(np.asarray(fill).flat[0])
                        sigma = np.nanmedian(np.where(valid, raw, np.nan), axis=0).T
                    observed_units.add(str(getattr(variable, "units", "")).strip())
                    data[band.name + "_wnum"] = np.asarray(x)
                    data[band.name + "_sigma"] = np.asarray(sigma)
                if len(observed_units) != 1:
                    raise ValueError(
                        f"{self.name} {self.platform}: noise units differ among bands: "
                        f"{sorted(observed_units)}."
                    )
                units = observed_units.pop()
        if source_platform != self.platform:
            raise ValueError(
                f"{self.name} {self.platform}: noise_file {source} identifies "
                f"platform {source_platform or 'unknown'}; provide a matching "
                "platform profile. NPZ profiles need a 'platform' field."
            )
        if self.noise_units is not None:
            units = self.noise_units
        if units.lower().replace(" ", "") not in (
            "mw/(m2srcm-1)", "mw/(m2*sr*cm-1)",
            "mw/(m2srcm^-1)", "mw",
        ):
            raise ValueError(
                f"{self.name} {self.platform}: noise radiance units {units!r} "
                "are unverified; supply an explicit noise_units='mW' only "
                "after confirming the profile is mW/(m2 sr cm-1)."
            )
        return data, source

    def add_noise(self, result: SRFM, rng: np.random.Generator | None = None) -> SRFM:
        if not isinstance(result, SRFM) or not hasattr(result, "uu"):
            raise ValueError(f"{self.name} noise needs a convolved SRFM radiance result.")
        if (getattr(result, "instrument_name", None) != self.name or
                getattr(result, "instrument_platform", None) != self.platform or
                getattr(result, "instrument_fov", None) != self.fov or
                getattr(result, "instrument_apodization", None) != self.apodization):
            raise ValueError(
                f"{self.name} {self.platform}: noise requires a channel result "
                "from the same instrument, platform and apodization."
            )
        if getattr(result, "noise_applied", False):
            raise ValueError(f"{self.name} {self.platform}: noise was already applied.")
        if self.values is None:
            raise ValueError("Noise requires instrument grid inputs.")
        expected = self.output_grid(self.values)
        if not np.array_equal(expected, np.asarray(result.wvnm)):
            raise ValueError(f"{self.name} {self.platform}: noise grid differs from selected native channels.")
        if rng is None:
            rng = np.random.default_rng(self.noise_seed)
        elif not isinstance(rng, np.random.Generator):
            raise TypeError("rng must be a numpy.random.Generator.")
        data, source = self._noise_arrays()
        sigma_by_band = {}
        for band in self.bands:
            x = np.asarray(data[band.name + "_wnum"], dtype=float)
            sigma = np.asarray(data[band.name + "_sigma"], dtype=float)
            if (x.shape != (band.count,) or sigma.shape != (self.fov_count, band.count)
                    or not np.allclose(x, band.channels, atol=1e-8, rtol=0)
                    or not np.all(np.isfinite(sigma) & (sigma >= 0))):
                raise ValueError(
                    f"{self.name} {self.platform} {band.name}: noise profile has "
                    "wrong channel coordinates, FOV count, or invalid amplitudes."
                )
            sigma_by_band[band.name] = sigma[self.fov - 1] * 1e-3
        radiance = np.asarray(result.uu)
        if radiance.shape[0] != len(expected) or not np.issubdtype(radiance.dtype, np.floating):
            raise ValueError("Radiance must be floating point and match the noise channel grid.")
        for spatial in np.ndindex(radiance.shape[1:]):
            draws = {}
            for band in self.bands:
                native = rng.normal(size=band.count) * sigma_by_band[band.name]
                if self.apodization == "Hamming":
                    native = np.convolve(native, HAMMING, mode="same")
                draws[band.name] = native
            radiance[(slice(None),) + spatial] += np.asarray([
                draws[self.bands[b].name][i] for b, i in self._selected
            ])
        if hasattr(result, "bbt"):
            result.calc_bbt()
        result.noise_applied = True
        result.noise_source = source
        result.processing_history = tuple(getattr(result, "processing_history", ())) + ({
            "operation": "add_noise", "instrument": self.name,
            "platform": self.platform, "fov": self.fov, "source": source,
            "apodization": self.apodization or "none",
        },)
        return result
