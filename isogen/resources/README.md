# Protein modification resources

`protein_modifications.json.gz` is a normalized runtime lookup generated from
the following local source files. Only fields needed for mass calculation are
retained.

- UniMod: `unimod2.obo.txt`, dated 2026-02-17, from
  <https://www.unimod.org/obo/unimod.obo>
- PSI-MOD: `PSI-MOD.obo`, data version 1.038.0 dated 2026-07-31, from
  <https://github.com/HUPO-PSI/psi-mod-CV>
- RESID: `RESIDUES.XML`, release 76.00 (2018-05-31), from
  <https://proteininformationresource.org/resid/>

The generated table contains 1,560 UniMod records, 1,639 PSI-MOD records, and
601 RESID records. Source names, accessions, mass deltas, formulas, residue
origins, and terminal/site specificity are retained where available. The
larger source databases, RESID images/models, and reference PDFs are not
distributed with IsoGen.

The table was generated on 2026-09-21. UniMod describes its database as public
domain under a copyleft permission requiring derivatives to retain the same
permissions. PSI-MOD is licensed under CC BY 4.0. RESID records retain their
source attribution and should be cited as described in its distributed README:
Garavelli, J.S., *Proteomics* 2004, 4, 1527-1533,
<https://doi.org/10.1002/pmic.200300777>. See the repository's
`THIRD_PARTY_NOTICES.md` for details.
