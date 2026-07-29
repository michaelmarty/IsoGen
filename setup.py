"""Setuptools hooks for bundling IsoGen's precompiled native library."""

import shutil
import sys
from pathlib import Path

from setuptools import Distribution, setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.build_py import build_py


ROOT = Path(__file__).resolve().parent


def platform_native_files():
    if sys.platform == "win32":
        return ("isogen.dll", "libmmd.dll", "svml_dispmd.dll")
    if sys.platform.startswith("linux"):
        return ("isogen.so",)
    raise RuntimeError("IsoGen wheels are currently supported only on Windows and Linux")


class BuildPyWithNativeLibrary(build_py):
    """Copy the current platform's native runtime into the wheel package."""

    def run(self):
        super().run()
        source_dir = ROOT / "bin"
        destination_dir = Path(self.build_lib) / "isogen" / "bin"
        destination_dir.mkdir(parents=True, exist_ok=True)

        for filename in platform_native_files():
            source = source_dir / filename
            if not source.is_file():
                raise FileNotFoundError("Required native library is missing: {}".format(source))
            shutil.copy2(source, destination_dir / filename)


class PlatformDistribution(Distribution):
    """Mark distributions as platform-specific despite having no Python extension."""

    def has_ext_modules(self):
        return True


class PlatformWheel(bdist_wheel):
    """Use a Python-agnostic ABI tag for the ctypes-based native library."""

    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        _, _, platform_tag = super().get_tag()
        return "py3", "none", platform_tag


setup(
    cmdclass={
        "build_py": BuildPyWithNativeLibrary,
        "bdist_wheel": PlatformWheel,
    },
    distclass=PlatformDistribution,
)
