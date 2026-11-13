"""Satellite-instrument processing strategies for SRFM results."""

from .base import ConvolvingInstrument, Instrument, NoiseInstrument
from .iasi import IASI, apply_banded_factor, normalize_iasi_instrument
from .registry import INSTRUMENT_NAMES, create_instrument, create_instrument_from_inputs

__all__ = [
    "ConvolvingInstrument",
    "IASI",
    "INSTRUMENT_NAMES",
    "Instrument",
    "NoiseInstrument",
    "apply_banded_factor",
    "create_instrument",
    "create_instrument_from_inputs",
    "normalize_iasi_instrument",
]
