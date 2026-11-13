"""Nominal FY-3D HIRAS-I full-resolution response and noise handling."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from .interferometer import Band, Interferometer


class HIRAS(Interferometer):
    name = "hiras"
    fov_count = 4
    # The supplied FY-3D DS_NEdN arrays have blank units in the granule and
    # conflicting product documentation.  The bundled amplitudes are usable
    # only with an explicit, scientifically checked hiras_noise_units='mW'.
    default_noise_resource = "hiras_fy3d_noise.npz"

    @classmethod
    def from_inputs(cls, values: Mapping[str, Any]) -> "HIRAS":
        return cls(
            platform="FY3D", fov=values.get("hiras_fov", 1),
            bands=(Band("LW", 648.75, 781), Band("MW1", 1208.75, 869),
                   Band("MW2", 2153.75, 637)),
            values=dict(values), apodization=values.get("apodization"),
            noise_file=values.get("hiras_noise_file"),
            noise_units=values.get("hiras_noise_units"),
            noise_seed=values.get("noise_seed"),
        )
