# Third-party notices

IsoGen's native libraries use third-party components. Before publishing a
binary distribution, the release manager is responsible for confirming that
the bundled binaries and notices correspond to the exact components used to
build that release.

## FFTW

The native IsoGen library is built with FFTW. The repository includes the FFTW
license and copyright notices in:

- `src/fftw/COPYING`
- `src/fftw/COPYRIGHT`

See <https://www.fftw.org/> for project and licensing information.

## Intel compiler runtime libraries

Windows wheels include `libmmd.dll` and `svml_dispmd.dll`, which are Intel
compiler runtime libraries required by `isogen.dll`. Intel documents these as
redistributable compiler libraries. Distribution is subject to the license
that accompanied the Intel compiler installation used to produce them.

Before publishing a Windows wheel, verify the applicable Intel license and
include its unmodified copyright notice and license terms. Intel's current
license index is available at:

<https://www.intel.com/content/www/us/en/developer/articles/license/end-user-license-agreement.html>
