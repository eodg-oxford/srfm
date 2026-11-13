"""Structural interfaces for independent instrument-processing operations."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

import numpy as np

if TYPE_CHECKING:
    from ..forward_model import SRFM


@runtime_checkable
class Instrument(Protocol):
    """Minimum identity shared by registered instrument implementations."""

    name: str


@runtime_checkable
class ConvolvingInstrument(Instrument, Protocol):
    """An instrument capable of convolving an SRFM radiance in place."""

    def output_grid(self, values: Mapping[str, Any]) -> np.ndarray:
        """Return the instrument channel grid selected by runner inputs."""

    def convolve(self, result: SRFM) -> SRFM:
        """Modify ``result`` in place and return the same instance."""


@runtime_checkable
class NoiseInstrument(Instrument, Protocol):
    """An instrument capable of adding noise to an SRFM result in place."""

    def add_noise(self, result: SRFM, **kwargs: Any) -> SRFM:
        """Modify ``result`` in place and return the same instance."""
