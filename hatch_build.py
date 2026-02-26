from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class SrfmBuildHook(BuildHookInterface):
    """Mark wheels as platform-specific and build the native extensions."""

    def dependencies(self) -> list[str]:
        """Ensure the build environment has the toolchain needed by f2py."""
        return ["numpy>=1.25", "meson>=1.3", "ninja>=1.11"]

    def initialize(self, version: str, build_data: dict[str, object]) -> None:
        build_data["pure_python"] = False
        build_data["infer_tag"] = True
        self._clean_artifacts()
        if self.target_name == "wheel":
            self._build_extensions()

    def clean(self, versions: list[str]) -> None:
        self._clean_artifacts()

    def _clean_artifacts(self) -> None:
        root = Path(self.root)
        sys.path.insert(0, str(root))
        try:
            from clean_srfm import clean
        finally:
            sys.path.pop(0)
        clean()

    def _build_extensions(self) -> None:
        """Invoke the shared build helper to compile all Fortran modules."""
        script = Path(self.root) / "build_extensions.py"
        if not script.exists():
            raise FileNotFoundError(f"Missing build helper: {script}")
        cmd = [sys.executable, str(script)]
        print(f"[build] running {' '.join(cmd)}")
        subprocess.run(cmd, cwd=self.root, check=True)
