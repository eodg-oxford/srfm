"""Nominal CrIS full-resolution channel response and J2 L1B noise."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from .interferometer import Band, Interferometer


class CrIS(Interferometer):
    name = "cris"
    fov_count = 9
    default_noise_resource = "cris_j2_noise.npz"

    @classmethod
    def from_inputs(cls, values: Mapping[str, Any]) -> "CrIS":
        platform = str(values.get("cris_platform", "J2")).upper().replace("-", "")
        if platform not in ("J2", "NOAA21"):
            raise ValueError(
                "CrIS native channel convolution currently supports J2/NOAA-21 "
                f"FSR; got cris_platform={platform!r}."
            )
        return cls(
            platform="J2", fov=values.get("cris_fov", 5),
            bands=(Band("LW", 648.75, 717), Band("MW", 1208.75, 869),
                   Band("SW", 2153.75, 637)),
            values=dict(values), apodization=values.get("apodization"),
            noise_file=values.get("cris_noise_file"),
            noise_units=values.get("cris_noise_units"),
            noise_seed=values.get("noise_seed"),
        )
