import numpy as np
import pytest

from srfm import disort_functions

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(("value", "expected"), [(2, 2), (3, 4), (16, 16)])
def test_update_maxcmu_enforces_even_stream_count(value, expected):
    """Verify stream-count updates produce the expected even value.

    DISORT requires an even number of computational streams.

    Args:
        value: Requested computational stream count.
        expected: Expected normalized stream count.
    """
    assert disort_functions.update_maxcmu(value) == expected


def test_update_dimension_helpers_follow_disort_flags():
    """Verify dimension helpers honor user-angle and flux-only flags.

    DISORT array sizes must track the active output configuration.
    """
    assert disort_functions.update_maxulv(1, False, 4) == 5
    assert disort_functions.update_maxulv(2, True, 4) == 2
    assert disort_functions.update_maxumu(1, False, 8) == 8
    assert disort_functions.update_maxumu(3, True, 8) == 3


def _valid_integrity_input():
    """Build a minimal valid mapping for DISORT integrity checks.

    Individual tests copy and alter this shared physical baseline.

    Returns:
        Complete small input mapping accepted by the integrity validator.
    """
    return dict(
        maxmom=4,
        maxcmu=4,
        maxumu=1,
        maxphi=1,
        ibcnd=0,
        onlyfl=False,
        dtauc=np.array([0.1, 0.2]),
        ssalb=np.array([0.5, 1.0]),
        temper=np.array([280, 270, 260]),
        wvnmlo=999,
        wvnmhi=1001,
        utau=np.array([0.3]),
        umu0=0.5,
        phi0=0,
        umu=np.array([0.8]),
        phi=np.array([0]),
        btemp=290,
        ttemp=250,
        temis=1,
    )


def test_disort_integrity_accepts_physically_valid_small_case():
    """Verify the integrity checker accepts a small physical case.

    The baseline mapping is reused by the invalid-field matrix below.
    """
    assert (
        disort_functions.test_disort_input_integrity(**_valid_integrity_input()) is True
    )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("maxmom", 2, "Maxmom"),
        ("maxcmu", 3, "even"),
        ("ssalb", np.array([1.1]), "ssalb"),
        ("wvnmlo", 1002, "Wvnmlo"),
        ("umu0", 2, "Umu0"),
        ("umu", np.array([0.8, 0.2]), "increasing"),
        ("phi", np.array([361]), "between"),
    ],
)
def test_disort_integrity_rejects_invalid_values(field, value, message):
    """Verify invalid DISORT fields report their specific constraint.

    Parameterization covers independent physical and dimensional failures.

    Args:
        field: Input field to replace.
        value: Invalid parameterized value.
        message: Expected diagnostic fragment.
    """
    values = _valid_integrity_input()
    values[field] = value
    with pytest.raises(ValueError, match=message):
        disort_functions.test_disort_input_integrity(**values)


def test_user_optical_depth_may_reach_total_column_depth():
    """Verify user optical depth may equal the total column depth.

    Equality represents the valid lower boundary rather than an overflow.
    """
    values = _valid_integrity_input()
    values["dtauc"] = np.array([0.1, 0.2])
    values["utau"] = np.array([0.3])
    assert disort_functions.test_disort_input_integrity(**values)
