"""Regression coverage for RFM-compatible output and scattering geometry."""

import numpy as np
import pytest

from srfm.rfm_levels import match_altitude, read_matching_grid, reconcile_levels

pytestmark = pytest.mark.unit


def test_collision_rebuilds_geometry_and_tracking_without_mutation():
    """Check canonical bounds and top-down tracking for the reported collision."""
    atmosphere = [0, 1, 2, 3, 4]
    requested = [0, 1, 2, 3, 4]
    boundaries = {"aerosol": (3.001, 3.101)}
    levels, tracking, effective = reconcile_levels(atmosphere, requested, boundaries)
    assert levels == [0, 1, 2, 3, float(np.float32(3.101)), 4]
    assert tracking == [None, "aerosol", None, None, None]
    assert effective["aerosol"] == (3, float(np.float32(3.101)))
    assert requested == atmosphere == [0, 1, 2, 3, 4]
    assert boundaries == {"aerosol": (3.001, 3.101)}


@pytest.mark.parametrize("altitude", [0, 1, 2])
def test_exact_matches(altitude):
    """Keep exactly matching atmospheric levels unchanged."""
    grid = np.array([0, 1, 2], dtype=np.float32)
    effective, updated = match_altitude(grid, altitude)
    assert effective == altitude
    np.testing.assert_array_equal(updated, grid)


@pytest.mark.parametrize("offset, snaps", [(-1, True), (0, True), (1, False)])
def test_single_precision_tolerance_boundary(offset, snaps):
    # An origin at zero makes the threshold itself representable in R4.
    """Check both adjacent R4 values and the inclusive tolerance threshold."""
    threshold = np.float32(0.001)
    altitude = threshold if offset == 0 else np.nextafter(
        threshold, np.float32(-np.inf if offset < 0 else np.inf)
    )
    effective, updated = match_altitude(np.array([0, 1], dtype=np.float32), altitude)
    assert (effective == 0) == snaps
    assert len(updated) == (2 if snaps else 3)


def test_decimal_collision_depends_on_single_precision_rounding():
    """Retain the single precision decision for the 3.001 km boundary."""
    grid = np.array([2, 3, 4], dtype=np.float32)
    assert float(np.float32(3.001) - np.float32(3)) < 0.001
    assert match_altitude(grid, 3.001)[0] == 3
    assert match_altitude(grid, np.nextafter(np.float32(3.001), np.float32(np.inf)))[0] > 3


def test_nearest_tie_uses_first_level_with_nonuniform_spacing():
    # The first tied level has a broad neighbour span and therefore snaps.
    """Select the first tied atmospheric level, as Fortran MINLOC does."""
    grid = np.array([0, 1000, 1001, 1002], dtype=np.float32)
    effective, updated = match_altitude(grid, 1000.5)
    assert effective == 1000
    np.testing.assert_array_equal(updated, grid)


@pytest.mark.parametrize("grid, altitude, expected", [
    ([0, 10, 11], 0.005, 0),
    ([0, 1, 11], 10.995, 11),
    ([0, 1, 11], 1.005, 1),
    ([0, 1, 11], 0.005, float(np.float32(0.005))),
])
def test_endpoint_and_interior_relative_spacing(grid, altitude, expected):
    """Use the local relative spacing at interior and endpoint levels."""
    effective, _ = match_altitude(np.array(grid, dtype=np.float32), altitude)
    assert effective == expected


def test_sequential_insertion_changes_later_tolerance():
    """Include earlier insertions when matching subsequent requests."""
    levels, tracking, _ = reconcile_levels([0, 1, 2], [0, 0.998, 1.0008, 2], {})
    assert levels == [0, float(np.float32(0.998)), float(np.float32(1.0008)), 2]
    assert tracking == [None] * 3
    # Without the earlier insertion, the second request snaps to 1 km.
    assert match_altitude(np.array([0, 1, 2], dtype=np.float32), 1.0008)[0] == 1


def test_valid_grid_and_multiple_requests_snapping_to_one_level():
    """Preserve valid levels and remove multiple requests for one effective level."""
    valid, _, _ = reconcile_levels([0, 1, 2], [0, 0.5, 1, 2], {})
    assert valid == [0, 0.5, 1, 2]
    levels, _, _ = reconcile_levels([0, 1, 2], [0, 0.0001, 0.0002, 1, 2], {})
    assert levels == [0, 1, 2]


def test_collapsed_layer_error_reports_original_and_effective_bounds():
    """Reject a collapsed physical layer with actionable geometry details."""
    with pytest.raises(ValueError, match="collapsed layer 'thin'.*original.*effective"):
        reconcile_levels([0, 10, 20], [0, 20], {"thin": (9.999, 10.001)})


@pytest.mark.parametrize("boundaries", [
    {"first": (0.5, 0.9995), "second": (1.0005, 1.5)},
    {"first": (0.5, 1), "second": (1, 1.5)},
    {"first": (0.5, 1.2), "second": (1, 1.5)},
])
def test_incompatible_layers_are_not_silently_discarded(boundaries):
    """Retain the runners' restriction against overlap and shared boundaries."""
    with pytest.raises(ValueError, match="first.*original.*effective.*second.*original.*effective"):
        reconcile_levels([0, 1, 2], [0, 2], boundaries)


@pytest.mark.parametrize("grid", [[0], [0, 0, 1], [1, 0], [0, np.nan]])
def test_invalid_atmospheric_grids(grid):
    """Reject grids that cannot be strictly increasing in RFM precision."""
    with pytest.raises(ValueError, match="atmospheric altitudes"):
        reconcile_levels(grid, [0, 1], {})


@pytest.mark.parametrize("subdivide", [False, True])
def test_first_atmospheric_grid_and_lay_preprocessing(
    tiny_altitude_grid, tmp_path, subdivide
):
    # The second file is intentionally absent: its composition grid must not
    # substitute for the first file's already established computational grid.
    """Use the computational grid and its LAY subdivision, not a composition grid."""
    driver = {"atmosphere": (tiny_altitude_grid.name, "composition.atm"),
              "flags": ("OPT", "NAD", "LEV") + (("LAY",) if subdivide else ())}
    grid = read_matching_grid(driver, {}, tmp_path)
    np.testing.assert_array_equal(grid, [0, .5, 1, 1.5, 2] if subdivide else [0, 1, 2])
    assert grid.dtype == np.float32


def test_jacobian_insertions_precede_output_levels(tiny_altitude_grid, tmp_path):
    """Replay Jacobian insertions before output-level matching."""
    driver = {"atmosphere": (tiny_altitude_grid,), "flags": ("NAD", "LEV", "JAC"),
              "jac": ("F11 0 .998 2",), "sfc": ("HGTSFC=1.5",)}
    grid = read_matching_grid(driver, {}, tmp_path)
    np.testing.assert_array_equal(grid, np.array([0, .998, 1, 2], dtype=np.float32))
    # Generated *SFC is after *LEV, so it cannot change these matching decisions.
    levels, _, _ = reconcile_levels(grid, [0, 1.0008, 2], {})
    assert levels[1] == float(np.float32(1.0008))


def test_existing_driver_order_and_contents_preserved(tiny_altitude_grid, tmp_path):
    """Respect existing section order and replace output levels only in memory."""
    from srfm.rfm_helper import _driver_with_output_levels

    path = tmp_path / "custom.drv"
    content = (f"*HDR\nexample\n*FLG\nOPT NAD LEV OBS\n*SPC\n999 1000 1\n"
               f"*GAS\nF11\n*ATM\n{tiny_altitude_grid}\n*TAN\n1\n"
               "*OBS\nHGTOBS=.998\n*LEV\n0 1.0008 2\n*END\n")
    path.write_text(content)
    grid = read_matching_grid({}, {"driver_path": path}, tmp_path)
    levels, _, _ = reconcile_levels(grid, [0, 1.0008, 2], {})
    assert levels[1] == float(np.float32(1.0008))
    records = _driver_with_output_levels(path, levels)
    assert records.index("*OBS") < records.index("*LEV")
    assert str(levels[1]) in records[records.index("*LEV") + 1]
    assert path.read_text() == content


@pytest.mark.parametrize("parameter", ["mass_loading", "n", "s_a_den", "v_den"])
def test_effective_thickness_precedes_concentration_conversion(
    tmp_path, parameter
):
    """Preserve each concentration parameterization using the effective thickness."""
    from srfm.layer import MieLayer
    from srfm.main import _reconcile_optical_layers
    from srfm.rfm_functions import write_atm_file

    grid = tmp_path / "grid.atm"
    write_atm_file({"HGT [km]": [2, 3, 4]}, grid)
    configured = MieLayer("aerosol", alt_low=3.001, alt_upp=3.101,
                          thick=.1, center_alt=3.051, rho=1800, r=.2, s=1.4,
                          **{parameter: 2.0})
    levels, tracking = _reconcile_optical_layers(
        [2, 3, 4], {"atmosphere": (grid,), "flags": ("NAD", "LEV")}, {}, tmp_path,
        {"aerosol": configured},
    )
    configured.nsv_or_ml()
    assert getattr(configured, parameter) == pytest.approx(2)
    assert configured.thick == levels[2] - levels[1]
    assert configured.center_alt == (levels[1] + levels[2]) / 2
    assert configured.mass_loading == pytest.approx(
        configured.size_distribution.volume_density * configured.thick * 1800e-6
    )
    assert tracking == [None, "aerosol", None]


def test_pressure_grid_lay_interpolates_altitude_in_log_pressure(tmp_path):
    """Use RFM log-pressure interpolation after pressure-grid subdivision."""
    path = tmp_path / "pressure.atm"
    path.write_text("3\n*PRE [mb]\n1000 500 250\n*HGT [km]\n0 5 10\n*END\n")
    grid = read_matching_grid({"atmosphere": (path,), "flags": ("LAY",)}, {}, tmp_path)
    np.testing.assert_allclose(grid, [0, 2.5, 5, 7.5, 10], atol=3e-6)
