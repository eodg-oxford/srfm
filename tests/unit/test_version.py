from importlib.metadata import version

import pytest

import srfm

pytestmark = pytest.mark.unit


def test_public_version_comes_from_distribution_metadata():
    """Verify the public version matches installed distribution metadata.

    This prevents package attributes from drifting away from release metadata.
    """
    assert srfm.__version__ == version("SRFM")
    assert srfm.version == srfm.__version__
