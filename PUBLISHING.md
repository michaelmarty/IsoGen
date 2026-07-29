# Building and publishing IsoGen

IsoGen contains precompiled native libraries, so Windows and Linux wheels must
be built and tested separately.

## Prepare a release

1. Update `isogen/_version.py`.
2. Update the changelog or release notes.
3. Confirm the native binaries in `bin/` were built from the intended commit.
4. Review `THIRD_PARTY_NOTICES.md` and include the exact licenses that apply to
   the FFTW and Intel runtime binaries being distributed.
5. Commit and tag the release, for example `v0.1.0`.

## Build locally

Install the release tools:

```shell
python -m pip install --upgrade build twine
```

Build a source distribution and a wheel for the current platform:

```shell
python -m build
python -m twine check dist/*
```

The expected wheel names are platform-specific:

- Windows: `isogen-<version>-py3-none-win_amd64.whl`
- Linux: initially `isogen-<version>-py3-none-linux_x86_64.whl`

For broad Linux compatibility, repair the Linux wheel in a compatible
manylinux environment. The build environment must provide `libfftw3.so.3` so
that `auditwheel` can bundle it:

```shell
python -m pip install auditwheel patchelf
auditwheel repair dist/isogen-*-linux_x86_64.whl --wheel-dir wheelhouse
```

## Test an artifact

Use a clean virtual environment and install the wheel itself, not the source
tree:

```shell
python -m venv wheel-test
wheel-test/Scripts/python -m pip install dist/isogen-*-win_amd64.whl
wheel-test/Scripts/python -c "import isogen; print(isogen.isodist(1000, isolen=8))"
wheel-test/Scripts/isogen dist PEPTIDE --isolen 8
```

On Linux, use `wheel-test/bin/python` and `wheel-test/bin/isogen`.

## Upload

TestPyPI is recommended for the first upload:

```shell
python -m twine upload --repository testpypi dist/*
```

After validating installation from TestPyPI:

```shell
python -m twine upload dist/*
```

The GitHub Actions workflow in `.github/workflows/publish.yml` builds both
platform wheels and publishes a GitHub release to PyPI through trusted
publishing. Configure the `pypi` environment and trusted publisher for this
repository on PyPI before publishing the first release.
