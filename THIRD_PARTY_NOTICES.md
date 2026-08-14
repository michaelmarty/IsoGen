# Third-party notices

IsoGen's native libraries use third-party components. 

## FFTW

The native IsoGen library is built with FFTW. The repository includes the FFTW
license and copyright notices in:

- `src/fftw/COPYING`
- `src/fftw/COPYRIGHT`

See <https://www.fftw.org/> for project and licensing information.

## Intel compiler runtime libraries

Windows wheels include `libmmd.dll` and `svml_dispmd.dll`, which are Intel
compiler runtime libraries required by `isogen.dll`. These are subject to the
Intel End User License Agreement (EULA) and redistribution terms. See:
<https://www.intel.com/content/www/us/en/developer/articles/license/end-user-license-agreement.html>

The PDF of the license is uploaded in `bin/`.
