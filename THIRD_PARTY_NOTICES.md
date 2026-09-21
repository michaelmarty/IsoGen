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

## Protein modification databases

`isogen/resources/protein_modifications.json.gz` contains a derived subset of
three protein-modification databases used for mass lookup:

- UniMod is described by its maintainers as a public-domain database under a
  copyleft permission allowing unrestricted redistribution and modification
  provided copies and derivatives retain the same permissions. See
  <https://www.unimod.org/unimod_help.html>.
- PSI-MOD is provided by the HUPO Proteomics Standards Initiative under the
  Creative Commons Attribution 4.0 International license. See
  <https://github.com/HUPO-PSI/psi-mod-CV>.
- RESID release 76.00 is attributed to John S. Garavelli and the Protein
  Information Resource. The distributed database README requests citation of
  Garavelli, J.S., *Proteomics* 2004, 4, 1527-1533,
  <https://doi.org/10.1002/pmic.200300777>.

IsoGen retains only names, accessions, mass deltas, formulas, and site/origin
metadata required for calculation; it does not redistribute the source
databases, images, models, or reference documents.
