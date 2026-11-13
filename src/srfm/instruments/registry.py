"""Explicit registry for built-in SRFM instrument implementations."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Callable

from .base import Instrument
from .iasi import IASI

InstrumentConstructor = Callable[..., Instrument]

_INSTRUMENTS: dict[str, InstrumentConstructor] = {
    "iasi": IASI,
}
INSTRUMENT_NAMES = frozenset(_INSTRUMENTS)


def _get_constructor(name: str) -> InstrumentConstructor:
    """Return the registered constructor for a case-insensitive name."""
    if not isinstance(name, str) or not name.strip():
        raise TypeError("Instrument name must be a non-empty string.")
    normalized_name = name.strip().lower()
    try:
        return _INSTRUMENTS[normalized_name]
    except KeyError as exc:
        available = ", ".join(sorted(_INSTRUMENTS))
        raise ValueError(
            f"Unknown instrument {name!r}. Available instruments: {available}."
        ) from exc


def create_instrument(name: str, **parameters: Any) -> Instrument:
    """Construct a registered instrument by its case-insensitive name."""
    return _get_constructor(name)(**parameters)


def create_instrument_from_inputs(
    values: Mapping[str, Any],
    selector: str = "convolve",
) -> Instrument | None:
    """Construct an instrument selected for one processing operation.

    Args:
        values: Validated runner input mapping.
        selector: Input field containing a registered instrument name or
            ``None``. The standard selectors are ``"convolve"`` and
            ``"add_noise"``.

    Returns:
        The configured instrument, or ``None`` when the operation is disabled.
    """
    name = values.get(selector)
    if name is None:
        return None

    constructor = _get_constructor(name)
    factory = getattr(constructor, "from_inputs", None)
    if factory is None:
        raise TypeError(
            f"Registered instrument {name!r} cannot be constructed from inputs."
        )
    return factory(values)
