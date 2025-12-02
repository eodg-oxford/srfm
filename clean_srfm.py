"""
Utility script to delete build artifacts produced by the SRFM toolchain.
Removes
- rfm outputs
- DISORT intermediates
- compiled modules (.so/.pyd)
- cache files
"""

from pathlib import Path
import shutil


REPO_ROOT = Path(__file__).resolve().parent


def iter_files(start_path: Path):
    for path in start_path.rglob("*"):
        if path.is_file():
            yield path


def should_remove_file(path: Path) -> bool:
    suffix_matches = {".asc", ".o", ".so", ".pyd", ".mod", ".pyc"}
    name_matches = {
        "rfm",
        "rfm.exe",
        "code.f",
        "combined.f90",
        "INTENSITY.dat",
        "rfm.log",
        "build.log",
    }
    return path.suffix in suffix_matches or path.name in name_matches


def clean() -> None:
    for file_path in iter_files(REPO_ROOT):
        if should_remove_file(file_path):
            file_path.unlink()

    removable_dirs = [
        REPO_ROOT / "src" / "srfm" / "__pycache__",
        REPO_ROOT / ".idea",
    ]
    for pycache in REPO_ROOT.rglob("__pycache__"):
        removable_dirs.append(pycache)

    for directory in removable_dirs:
        if directory.is_dir():
            shutil.rmtree(directory)

    print("Successfully cleaned srfm.")


if __name__ == "__main__":
    clean()
