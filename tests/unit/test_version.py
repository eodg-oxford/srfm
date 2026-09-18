from importlib.metadata import version
from pathlib import Path
import tomllib

import pytest

import srfm

pytestmark = pytest.mark.unit


def test_public_version_comes_from_distribution_metadata():
    """Verify the public version matches installed distribution metadata.

    This prevents package attributes from drifting away from release metadata.
    """
    assert srfm.__version__ == version("SRFM")
    assert srfm.version == srfm.__version__


def test_release_metadata_and_license_file_use_gpl_v3():
    """Verify distributions declare and include the project-wide GPL terms."""
    repository = Path(__file__).resolve().parents[2]
    with (repository / "pyproject.toml").open("rb") as handle:
        metadata = tomllib.load(handle)

    assert metadata["project"]["license"] == "GPL-3.0-only"
    assert "LICENSE" in metadata["project"]["license-files"]
    license_text = (repository / "LICENSE").read_text(encoding="utf-8")
    assert "GNU GENERAL PUBLIC LICENSE" in license_text
    assert "Version 3, 29 June 2007" in license_text
