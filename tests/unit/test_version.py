from importlib.metadata import version

import pytest

import srfm


pytestmark = pytest.mark.unit


def test_public_version_comes_from_distribution_metadata():
    assert srfm.__version__ == version("SRFM")
    assert srfm.version == srfm.__version__
