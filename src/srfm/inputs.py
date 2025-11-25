"""This module defines the srfm input class and it's methods.

- Name: inputs
- Parent package: srfm
- Author: Antonin Knizek
- Contributors:
- Date: 25 Nov 2025
"""

from __future__ import annotations
from pathlib import Path
import importlib.util
import uuid
from mergedeep import merge

class Inputs:
    """Class that collects inputs for the main srfm run.

    Has various methods to read and combine inputs from different sources.

    """

    def __init__(self, **kwargs):
        """Init function, set one explicit parameter values, which is a dict containing
        the actual parameters for srfm.

        """
        self.values = {}
        for key, val in kwargs.items():
            setattr(self, key, val)

    def read_srfm_drv(self, drv: str | Path) -> None:
        """Read srfm driver table and load inputs into the values dictionary.

        Args:
            drv (str, Path): Path to driver table.

        Raises:
            TypeError: Raised when input path format not recognized.
            FileNotFoundError: Raised when driver table not found.

        """
        if isinstance(drv, str):
            drv = Path(drv).expanduser().resolve()
        elif isinstance(drv, Path):
            pass
        else:
            raise TypeError("drv must be str or pathlib.Path.")

        if not drv.is_file():
            raise FileNotFoundError(f"Driver path not found: {drv}.")

        module_name = f"_driver_table_{uuid.uuid4().hex}"
        spec = importlib.util.spec_from_file_location(module_name, drv)
        if spec is None or spec.loader is None:
            raise ImportError(f"Unable to load driver table from {filepath}")

        module = importlib.util.module_from_spec(spec)
        loader = spec.loader
        loader.exec_module(module)

        self.values = getattr(module, "inputs")
        return

    def read_oxharp_drv(self, drv: str | Path) -> None:
        """Read oxharp driver table and load inputs into the values dictionary.

        Args:
            drv (str, Path): Path to driver table.

        Raises:
            TypeError: Raised when input path format not recognized.
            FileNotFoundError: Raised when driver table not found.

        """

        if isinstance(drv, str):
            drv = Path(drv).expanduser().resolve()
        elif isinstance(drv, Path):
            pass
        else:
            raise TypeError("drv must be str or pathlib.Path.")

        if not drv.is_file():
            raise FileNotFoundError(f"Driver path not found: {drv}.")

        module_name = f"_driver_table_{uuid.uuid4().hex}"
        spec = importlib.util.spec_from_file_location(module_name, drv)
        if spec is None or spec.loader is None:
            raise ImportError(f"Unable to load driver table from {filepath}")

        module = importlib.util.module_from_spec(spec)
        loader = spec.loader
        loader.exec_module(module)

        _x = getattr(module, "STATE")
        _b = getattr(module, "ANCILLARY")
        
        self.values = merge({}, _x, _b)
        return
