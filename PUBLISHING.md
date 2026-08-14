# Building and publishing IsoGen

IsoGen uses scikit-build-core to compile its native C library while each wheel
is built. Windows, Linux, macOS Intel, and macOS Apple Silicon wheels must be
built and tested on their target platforms.

## Prepare a release

1. Update `isogen/_version.py`.
2. Update the changelog or release notes.
3. Confirm the native sources and bundled model files are from the intended commit.
4. Review `THIRD_PARTY_NOTICES.md` and include the exact licenses that apply to
   the FFTW and Intel runtime binaries being distributed.
5. Commit the release changes and push them to GitHub.

## Publish with GitHub Actions

Run the **Build and publish** workflow from the Actions tab on the commit you
want to release. The workflow reads the version from `isogen/_version.py`,
builds and tests all platform wheels plus the source distribution, creates the
matching `v<version>` tag, and publishes a GitHub release containing those
artifacts. It fails without changing an existing release if that tag has
already been released.

Set **Also publish the release artifacts to PyPI?** to `true` to publish the
same artifacts through PyPI trusted publishing after the GitHub release is
created. Configure the `pypi` environment and trusted publisher for this
repository on PyPI before using that option.

## Build locally

Linux CMake builds enable `ISOGEN_STATIC_GNU_RUNTIME` by default. This passes
`-static-libstdc++ -static-libgcc` when linking the native artifacts so a
binary built with a recent GCC does not require that GCC version's
`GLIBCXX_*` symbols on the target system. Distribution maintainers may opt
out with `-DISOGEN_STATIC_GNU_RUNTIME=OFF` when the runtime dependency is
managed by the distribution.

Install the release tools:

```shell
python -m pip install --upgrade build twine
```

Delete all previous build output before creating a release. This prevents an
older platform wheel from being uploaded with the new release:

```shell
python -c "import shutil; [shutil.rmtree(path, ignore_errors=True) for path in ('build', 'dist', 'wheelhouse')]"
```

Build a source distribution and a wheel for the current platform. Wheel builds
run CMake automatically and require a compiler plus the FFTW development
library:

```shell
python -m build
python -m twine check dist/*
```

The expected wheel names are platform-specific:

- Windows: `pyisogen-<version>-py3-none-win_amd64.whl`
- Linux: initially `pyisogen-<version>-py3-none-linux_x86_64.whl`
- macOS Intel: `pyisogen-<version>-py3-none-macosx_*_x86_64.whl`
- macOS Apple Silicon: `pyisogen-<version>-py3-none-macosx_*_arm64.whl`

For broad Linux compatibility, repair the Linux wheel in a compatible
manylinux environment. The build environment must provide `libfftw3.so.3` so
that `auditwheel` can bundle it:

```shell
python -m pip install auditwheel patchelf
auditwheel repair dist/pyisogen-*-linux_x86_64.whl --wheel-dir wheelhouse
```

Before committing a rebuilt Linux library, confirm that it has no dynamic GNU
runtime dependency:

```shell
readelf --dynamic bin/isogen.so | grep -E 'libstdc\+\+|libgcc_s'
readelf --version-info bin/isogen.so | grep GLIBCXX_
```

Both commands should produce no output.

On macOS, install FFTW before building and use `delocate-wheel` to bundle and
relink its runtime library:

```shell
brew install fftw
python -m pip install delocate
mkdir wheelhouse
delocate-wheel --require-archs "$(uname -m)" -w wheelhouse -v dist/*.whl
```

## Test an artifact

Use a clean virtual environment and install the wheel itself, not the source
tree:

```shell
python -m venv wheel-test
wheel-test/Scripts/python -m pip install dist/pyisogen-*-win_amd64.whl
wheel-test/Scripts/python -c "import isogen; print(isogen.isodist(1000, isolen=8))"
wheel-test/Scripts/python -c "import isogen; print(isogen.isodist('C6H12O6', type='ATOM', isolen=8))"
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

The GitHub Actions workflow in `.github/workflows/publish.yml` builds Windows,
Linux, macOS Intel, and macOS Apple Silicon wheels. A manually dispatched run
creates the GitHub release and can optionally publish the same files to PyPI.
