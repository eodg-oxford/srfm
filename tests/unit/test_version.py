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


def test_release_metadata_is_version_1_2_0():
    """Verify the source metadata identifies the 1.2.0 release."""
    project_file = Path(__file__).resolve().parents[2] / "pyproject.toml"
    with project_file.open("rb") as handle:
        metadata = tomllib.load(handle)

    assert metadata["project"]["version"] == "1.2.0"
