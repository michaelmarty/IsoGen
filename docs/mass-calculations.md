# Mass calculations

IsoGen exposes average and monoisotopic neutral-mass functions for proteins,
RNA, DNA, and elemental formulas.

## Elemental formulas

Calculate a formula's light-isotope mass or matching mass axis:

```python
formula_mass = isogen.calc_atom_monoisotopic_mass("C6H12O6")
formula_axis = isogen.calc_atom_mass_axis("C6H12O6", isolen=8)
```

The formula parser supports the same 109 elements as the native ATOM FFT
calculation. The count vector is available from
`isogen.isogenwrapper.atom_formula_to_vector`.

## Proteins

An intact protein uses `H2O` terminal composition by default:

```python
average = isogen.calc_pep_mass("PEPTIDE")
monoisotopic = isogen.calc_pep_monoisotopic_mass("PEPTIDE")
```

### Modified proteins (ProForma)

Protein mass functions accept the mass-bearing parts of the current
[HUPO-PSI ProForma notation](https://github.com/HUPO-PSI/ProForma). UniMod,
PSI-MOD, and RESID names and accessions are resolved from databases bundled
with IsoGen. Numeric shifts, elemental formulas, terminal modifications,
labile or unlocalized modifications, and global fixed modifications are also
supported:

```python
oxidized = isogen.calc_pep_monoisotopic_mass("EM[Oxidation]E")
psi_mod_oxidized = isogen.calc_pep_monoisotopic_mass(
    "EM[MOD:00719]E"
)
resid_oxidized = isogen.calc_pep_monoisotopic_mass(
    "EM[RESID:AA0581]E"
)
n_terminal = isogen.calc_pep_monoisotopic_mass("[Acetyl]-PEPTIDE")
fixed_cysteine = isogen.calc_pep_monoisotopic_mass(
    "<[Carbamidomethyl]@C>ACDC"
)
```

`J`, `O`, and `U` have defined residue masses. `B` and `Z` use the arithmetic
midpoint of their D/N and E/Q possibilities. Because `X` has no unique mass,
it must carry a complete signed mass gap, for example
`RTAAX[+367.0537]WT`; bare `X` raises `ValueError`.

XL-MOD, GNO/GNOme, cross-links, branched peptides, and multi-peptidoform
expressions are not yet supported by the strict mass functions. Global isotope
replacement, such as `<13C>PEPTIDE`, and explicit isotope formulas are supported
when their composition and isotope masses can be resolved.

`isodist` includes resolvable modification composition by default:

```python
modified = isogen.isodist("EM[Oxidation]E")
legacy = isogen.isodist("EM[Oxidation]E", use_modifications=False)
labeled = isogen.isodist("K[UNIMOD:259]")
formula_label = isogen.isodist("A[Formula:H-2[2H2]]")
```

`use_modifications=False` restores the legacy intensity calculation. Known
modification mass shifts still affect the axis in both modes. A numeric mass
shift alone cannot determine composition: prediction prints a warning to
standard error and ignores it for intensities. Unknown or unsupported descriptors
are skipped individually, allowing supported modifications on the same peptide
to remain active. The standalone mass functions retain their strict validation.

Unimod, PSI-MOD, and site-specific RESID formulas are resolved offline. Unimod
molecular subunits are expanded into atoms; common `Glycan:` compositions are
also supported by distribution prediction. Fixed isotope atoms contribute their
exact mass but no natural-abundance envelope width. Their precomputed FFT lives
in a separate native object, leaving the original natural-abundance tables intact.
Partial isotope enrichment is not supported.

Composition-aware distributions use light-isotope mass origins for elements such
as selenium and iron, matching the native general-element FFT. These can differ
from CV monoisotopic masses defined using the most abundant isotope. The legacy
residue/CV mass APIs retain their existing conventions.

For composition-changing modifications, bundled NN methods warn and use FFT.
BRAIN handles CHNOS plus phosphorus, and warns before using FFT for other
elements or fixed isotope labels. Custom NN models warn and ignore composition
changes in intensities because their residue-count input cannot represent them.

Use the dedicated helpers when parsing is useful independently:

```python
mass = isogen.calc_proforma_mass("EM[UNIMOD:35]E")
sequence = isogen.strip_proforma("EM[UNIMOD:35]E")  # "EME"
```

The supported `ion_type` values describe neutral terminal compositions:

| Ion type | Sequence to supply | Terminal shift from residue sum |
| --- | --- | --- |
| `H2O` | Full protein | +H2O |
| `a` | N-terminal fragment | -CO |
| `a+1` | N-terminal fragment | -CO+H |
| `b` | N-terminal fragment | none |
| `c` | N-terminal fragment | +NH3 |
| `x` | C-terminal fragment | +CO2 |
| `x+1` | C-terminal fragment | +CO2+H |
| `y` | C-terminal fragment | +H2O |
| `y-1` | C-terminal fragment | +H2O-H |
| `z` | C-terminal fragment | +H2O-NH3 (Pyteomics `z`) |
| `z'` | C-terminal fragment | +H2O-NH2 (one H above `z`) |

`z+1`, `z•`, `z·`, and `z.` are accepted as aliases for `z'`.

`calc_pep_fragments` can select ion series by fragmentation method:

| `fragmentation_type` | Ion series |
| --- | --- |
| `CID`, `HCD`, `SID`, `IRMPD` | `b`, `y` |
| `ETD`, `ECD` | `c`, `z'` |
| `EThcD`, `BYCZ*` | `b`, `y`, `c`, `z'` |
| `UVPD` | `a`, `b`, `c`, `x`, `y`, `z'` |
| `UVPD4` | `a`, `a+1`, `x+1`, `y-1` |
| `UVPD6` | `a`, `a+1`, `x+1`, `x`, `y-1`, `z'` |
| `UVPD9` | `a`, `a+1`, `b`, `c`, `x`, `x+1`, `y`, `y-1`, `z'` |

```python
b6_mass = isogen.calc_pep_monoisotopic_mass("PEPTID", ion_type="b")
y6_mass = isogen.calc_pep_monoisotopic_mass("EPTIDE", ion_type="y")
fragments = isogen.calc_pep_fragments("PEPTIDE")
# {"b1": ..., "b2": ..., ..., "y1": ..., "y2": ..., ...}
ecd_fragments = isogen.calc_pep_fragments("PEPTIDE", ion_types=("c", "z'"))
uvpd_fragments = isogen.calc_pep_fragments("PEPTIDE", fragmentation_type="UVPD9")

modified = isogen.calc_pep_fragments("EM[Oxidation]E")
ambiguous = isogen.calc_pep_fragments(
    "AS[#g1]T[Phospho#g1]K", ambiguous_rule="both"
)
# A cleavage separating the candidate sites produces "b2#1" and "b2#2".
```

Supplying `ion_types` explicitly overrides `fragmentation_type`.

Localized, terminal, global fixed, region-localized, labile, and unlocalized
ProForma modifications are supported. The default `ambiguous_rule="reject"`
omits any fragment with more than one possible mass. With
`ambiguous_rule="both"`, each distinct possibility is returned as a float
under a numbered key such as `b2#1`, ordered by increasing mass. Ambiguous
`B` and `Z` residues follow the same rule. `J` remains unambiguous by mass
because isoleucine and leucine have equal residue masses.

These are neutral masses. Charge and proton/adduct masses are not applied.
When used through `isodist`, `ion_type` changes the mass-axis origin. The
sequence intensity calculation retains its standard intact-sequence terminal
composition.

### Fragment isotope envelopes

`calc_pep_fragment_isodists` generates the selected peptide fragment masses
and isotope envelopes together:

```python
batch = isogen.calc_pep_fragment_isodists(
    "S[Acetylation]HHS", fragmentation_type="ETD", isolen=64
)
for label, mass, intensities in zip(batch.labels, batch.masses, batch.intensities):
    print(label, mass, intensities)
```

`batch.labels`, `batch.masses`, and rows of `batch.intensities` are aligned.
Masses are neutral monoisotopic origins and intensities have shape
`(number_of_fragments, isolen)`, normalized per fragment to its base peak.
The function supports the same terminal-ion selection as
`calc_pep_fragments`, but requires `monoisotopic=True` and an `isolen` from 1
through 128.

Exact fragment envelopes require an unambiguous natural-abundance CHNOS
composition. A supported ProForma modification is applied only to the
fragments that retain it. Mass-only, fixed-isotope, unsupported-element,
unlocalized, and ambiguously localized modifications raise `ValueError`
rather than returning an incomplete envelope.

## RNA and DNA termini

Nucleic-acid calculations default to a 3'-OH and a 5'-monophosphate:

```python
rna_mass = isogen.calc_rna_monoisotopic_mass("AUGC")
dna_mass = isogen.calc_dna_monoisotopic_mass("ATGC")
```

`threeend` accepts `OH` or an unrecognized/no-adjustment value. `fiveend`
accepts:

| Value | Meaning |
| --- | --- |
| `OH` | 5'-hydroxyl |
| `MP` | 5'-monophosphate (default) |
| `TP` | 5'-triphosphate |

```python
triphosphate = isogen.calc_rna_monoisotopic_mass(
    "AUGC",
    threeend="OH",
    fiveend="TP",
)
```

RNA `a`, `b`, `c`, `d`, `w`, `x`, `y`, and `z` fragment-series names are not
accepted by `ion_type`; that keyword applies only to peptides. Pass a manually
truncated RNA sequence and select from the terminal chemistries above, or pass
a known neutral fragment mass as numeric RNA input. Numeric input uses the RNA
averagine intensity model.

When `threeend` or `fiveend` is passed through `isodist`, it adjusts the mass
axis but does not alter the sequence-model intensity vector's terminal
composition.

## Standalone mass axes

Create an axis from a known first mass:

```python
axis = isogen.calc_mass_axis(1000.0, isolen=8, isotope_spacing=1.0033)
```

Or calculate the origin from a sequence:

```python
peptide_axis = isogen.calc_pep_mass_axis("PEPTIDE", isolen=8)
rna_axis = isogen.calc_rna_mass_axis("AUGC", isolen=8)
dna_axis = isogen.calc_dna_mass_axis("ATGC", isolen=8)
formula_axis = isogen.calc_atom_mass_axis("C6H12O6", isolen=8)
```

`gen_mass_axis` dispatches between numeric, sequence, and formula input in the
same way as `isodist`.
