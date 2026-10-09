"""Regression coverage for the October 2026 bundled ARIA replacement."""

from hashlib import sha256
from importlib.resources import files
import json
from pathlib import Path
import re

import numpy as np
import pytest

from srfm.ARIA_module import RI, ReadError, get_ri_filepathname, read_ri_file
from srfm.optical_properties import get_ri

pytestmark = pytest.mark.unit

LEGACY_MANIFEST = json.loads(
    (Path(__file__).parents[1] / "fixtures/unit/aria_legacy_manifest.json").read_text(
        encoding="utf-8"
    )
)
BUNDLED_ARIA_DIRECTORY = Path(str(files("srfm.data") / "ARIA"))
BUNDLED_RI_FILES = sorted(BUNDLED_ARIA_DIRECTORY.rglob("*.ri"))

# Each tuple records wavelength, column, old value, and supplied replacement.
# These are the only approved numerical changes among the 450 original files.
PETERSON_QUARTZ_CHANGES = {
    "quartz_E_Peterson_1969.ri": [(17.8, "k", 0.07604, 0.07904)],
    "quartz_O_Peterson_1969.ri": [
        (9.0, "n", 2.59701, 0.17463),
        (9.4, "n", 1.37131, 6.38644),
        (9.4, "k", 4.87845, 1.37131),
        (22.3, "k", 3.03182, 3.03132),
        (25.4, "n", 5.99961, 5.99862),
    ],
}


def test_bundled_database_has_unique_filenames():
    """Ensure basename lookup cannot silently select another dataset."""
    assert len(BUNDLED_RI_FILES) == 677
    assert len({path.name for path in BUNDLED_RI_FILES}) == len(BUNDLED_RI_FILES)
    assert len(LEGACY_MANIFEST) == 450


@pytest.mark.parametrize(
    "reference", LEGACY_MANIFEST, ids=lambda item: Path(item["old_path"]).name
)
def test_legacy_filename_preserves_original_numeric_data(reference):
    """Compare every old dataset against its verified replacement.

    Fingerprints include every stored numeric column and row, including files
    that contain only n or k and cannot supply a complete Mie refractive index.
    For the two Peterson quartz files, assert each approved change, then undo
    those changes in memory to compare all remaining values with the old file.

    Args:
        reference: Numeric fingerprint captured from the old database.
    """
    legacy_filename = Path(reference["old_path"]).name
    replacement = Path(get_ri_filepathname(legacy_filename))
    assert replacement.name == reference["new_filename"]
    assert replacement == Path(get_ri_filepathname(reference["new_filename"]))
    content = replacement.read_text(encoding="utf-8")
    columns = re.search(
        r"^#\s*FORMAT\s*=\s*(.*)$", content, re.MULTILINE | re.IGNORECASE
    ).group(1).lower().split()
    assert columns == reference["columns"]
    numeric_data = np.loadtxt(replacement, dtype="<f8", ndmin=2, encoding="utf-8")
    assert list(numeric_data.shape) == reference["shape"]
    for wavelength, column, old_value, new_value in PETERSON_QUARTZ_CHANGES.get(
        legacy_filename, []
    ):
        row_indices = np.flatnonzero(
            numeric_data[:, columns.index("wavl")] == wavelength
        )
        assert len(row_indices) == 1
        row_index = row_indices[0]
        column_index = columns.index(column)
        assert numeric_data[row_index, column_index] == new_value
        numeric_data[row_index, column_index] = old_value
    assert sha256(numeric_data.tobytes()).hexdigest() == reference["numeric_sha256"]

    if (
        "selection_sha256" in reference
        and legacy_filename not in PETERSON_QUARTZ_CHANGES
    ):
        refractive_indices = RI()
        refractive_indices.read(replacement)
        for mode, column in [("wavelength", "wavl"), ("wavenumber", "wavn")]:
            coordinates = refractive_indices.data[column]
            target_grid = np.linspace(min(coordinates), max(coordinates), 7)
            interpolated = np.asarray(
                refractive_indices.select(target_grid, mode=mode), dtype="<f8"
            )
            assert (
                sha256(interpolated.tobytes()).hexdigest()
                == reference["selection_sha256"][mode]
            )


@pytest.mark.parametrize("refractive_index_path", BUNDLED_RI_FILES, ids=lambda path: path.name)
def test_current_database_files_resolve_and_preserve_available_columns(
    refractive_index_path,
):
    """Check current names, new datasets, and the incomplete-data diagnostic.

    Args:
        refractive_index_path: File in the supplied replacement database.
    """
    assert (
        Path(get_ri_filepathname(refractive_index_path.name))
        == refractive_index_path
    )
    columns = re.search(
        r"^#\s*FORMAT\s*=\s*(.*)$",
        refractive_index_path.read_text(encoding="utf-8"),
        re.MULTILINE | re.IGNORECASE,
    ).group(1).lower().split()
    refractive_indices = RI()
    if not {"n", "k"}.issubset(columns):
        with pytest.raises(ReadError, match="FORMAT requires n and k"):
            refractive_indices.read(refractive_index_path)
        return
    refractive_indices.read(refractive_index_path)
    numeric_data = np.loadtxt(refractive_index_path, ndmin=2, encoding="utf-8")
    for column_index, column in enumerate(columns):
        np.testing.assert_array_equal(
            refractive_indices.data[column], numeric_data[:, column_index]
        )


@pytest.mark.parametrize(
    "legacy_filename",
    [
        "eyjafjallajokull-ash_Reed.ri",
        "ICE_Warren_2008.ri",
        "H2SO4_75_Palmer_1975.ri",
        "H2O_263K_Rowe_2020.ri",
        "quartz100_Henning_1997.ri",
        "malic-acid_Laskina_2014.ri",
    ],
)
def test_original_files_and_replacements_produce_identical_mie_indices(
    legacy_filename, legacy_aria_directory
):
    """Compare real old files and replacements through the Mie index loader.

    Args:
        legacy_filename: Original basename retained in the regression archive.
        legacy_aria_directory: Directory containing six original ARIA files.
    """
    reference = next(
        entry for entry in LEGACY_MANIFEST
        if Path(entry["old_path"]).name == legacy_filename
    )
    old_path = legacy_aria_directory / legacy_filename
    assert sha256(old_path.read_bytes()).hexdigest() == reference["old_file_sha256"]
    wavelengths = 10000.0 / np.array([999.5, 1000.0, 1000.5])
    old_real, old_extinction = read_ri_file(old_path, wave=wavelengths)
    for composition in (legacy_filename, reference["new_filename"]):
        for mode in ("wavelength", "wavenumber"):
            original_arrays = read_ri_file(old_path, mode=mode)
            replacement_arrays = RI().load_refractive_indices(composition, mode=mode)
            for original, replacement in zip(original_arrays, replacement_arrays):
                np.testing.assert_array_equal(replacement, original)
        replacement_indices = get_ri(composition, wave=wavelengths, wave_size=3)
        np.testing.assert_array_equal(
            replacement_indices, old_real - 1j * old_extinction
        )
