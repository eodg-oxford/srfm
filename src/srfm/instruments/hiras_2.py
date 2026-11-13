"""FY-3E/F/H HIRAS-II native full-resolution responses and noise."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from .interferometer import Band, Interferometer


class HIRAS2(Interferometer):
    name = "hiras-2"
    fov_count = 9

    @classmethod
    def from_inputs(cls, values: Mapping[str, Any]) -> "HIRAS2":
        platform = str(values.get("hiras2_platform", "FY3E")).upper().replace("-", "")
        if platform not in ("FY3E", "FY3F", "FY3H"):
            raise ValueError(
                "HIRAS-2 platform must be FY3E, FY3F or FY3H; "
                f"got hiras2_platform={platform!r}."
            )
        if platform == "FY3H":
            bands = (
                Band("LW", 648.75, 821, "table"),
                Band("MW1", 1159.375, 1292, "table"),
                Band("MW2", 1964.375, 940, "table"),
            )
        else:
            response = "table" if platform == "FY3F" else "sinc"
            bands = (
                Band("LW", 648.75, 834, "table"),
                Band("MW1", 1167.5, 1207, response),
                Band("MW2", 1919.375, 1012, response),
            )
        instrument = cls(
            platform=platform, fov=values.get("hiras2_fov", 1), bands=bands,
            values=dict(values), apodization=values.get("apodization"),
            noise_file=values.get("hiras2_noise_file"),
            noise_units=values.get("hiras2_noise_units"),
            noise_seed=values.get("noise_seed"),
        )
        instrument.default_noise_resource = (
            "hiras2_fy3e_noise.npz" if platform == "FY3E" else ""
        )
        return instrument
