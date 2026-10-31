# Isotope distributions

The top-level `isodist` function combines an intensity vector from the native
library with a matching neutral-mass axis.

```python
isogen.isodist(
    input,
    type="PEPTIDE",
    isolen=128,
    method="FFT",
    **mass_kwargs,
)
```

For peptide, RNA, and DNA inputs, `method` accepts `FFT`, `BRAIN`, or `NN`.
FFT and BRAIN are absolute calculations; BRAIN uses a polynomial recurrence,
while NN uses a bundled trained model. Method names are case-insensitive in
the Python API.

## Numeric mass input

When `input` is numeric, it becomes the first point on the mass axis:

```python
distribution = isogen.isodist(25_000, isolen=32)
assert distribution[0, 0] == 25_000
```

Numeric input selects an averagine-style intensity model. It does not infer an
elemental composition from a sequence.

## Sequence input

String input selects the corresponding sequence model. The mass-axis origin is
calculated from the sequence:

```python
distribution = isogen.isodist("PEPTIDE", type="PEPTIDE")
assert distribution[0, 0] == isogen.calc_pep_monoisotopic_mass("PEPTIDE")
```

RNA accepts `T` as uracil, and DNA accepts `U` as thymine for mass
calculations.

## Elemental-formula input

Select `type="ATOM"` to calculate a direct isotope distribution from a
formula:

```python
distribution = isogen.isodist(
    "C6H12O6",
    type="ATOM",
    isolen=32,
)
```

ATOM uses the FFT method only. Its first mass is the formula's light-isotope
mass, and its intensity vector is base-peak normalized. Formulas use standard
element symbols followed by optional nonnegative integer counts, such as
`H2O`, `C6H12O6`, or `NaCl`.

### Partially labeled (enriched) samples

Pass `heavy_isotopes` to describe a sample enriched in a heavy isotope:

```python
distribution = isogen.isodist(
    "C6H12O6",
    type="ATOM",
    isolen=32,
    heavy_isotopes={"2H": 0.05},  # 5% of the hydrogens are deuterium
)
```

Keys are a mass number and element symbol (`"2H"` or its alias `"D"`,
`"13C"`, `"15N"`, `"18O"`, or any tabulated isotope). Each value is the
fraction of that element's atoms that are the heavy isotope. The element's
single-atom distribution becomes the linear combination
`(1 - f) * natural + f * heavy`, where the heavy isotope is a delta function;
that mixed distribution is Fourier transformed and raised to the atom count in
place of the precomputed natural-abundance transform. Several elements can be
labeled at once, and several isotopes of one element (such as `17O` and `18O`)
are allowed when their fractions sum to at most 1.

Fractions must be between 0 and 1, and each labeled element must appear in the
formula; otherwise `ValueError` is raised. `heavy_isotopes` applies only to
`ATOM` formula inputs. The mass axis is unchanged and still starts at the
light-isotope mass. The natural abundance of the heavy isotope is part of the
`(1 - f)` natural fraction, so `f` is not added on top of it. Enriched elements
bypass the precomputed transform, so they cost slightly more per call.

### Measuring enrichment

`calc_heavy_percent` does the reverse: given a measured distribution, an
elemental formula, and the heavy isotope to test, it returns the percent
enrichment and the fitted distribution.

```python
measured = isogen.isodist("C6H12O6", type="ATOM", isolen=24, dist_only=True,
                          heavy_isotopes={"2H": 0.05})
percent, fit = isogen.calc_heavy_percent(measured, "C6H12O6", "2H")  # ~5.0
```

`dist` is a 1-D intensity vector (or an `(m, 2)` array whose second column is
used) that starts at the light-isotope peak and has one value per nominal mass
unit; its length is the fitted window. The overall intensity scale is
arbitrary. Each peak is a degree-n polynomial in the enrichment `f` (n is the
atom count of the element), because the element contributes the n-fold
convolution of `(1 - f) * natural + f * heavy` with the fixed distribution of
all other elements. The scale and `f` are fitted by bounded least squares
(`scipy.optimize.least_squares` with an analytic Jacobian), started from the
best point of a coarse scan that includes the closed-form mean-shift estimate.
The standard error comes from the Gauss-Newton covariance. `fit_dist` is
`isodist(formula, heavy_isotopes={heavy_atom: f})` at the fitted `f`. Use
`weights="poisson"` to weight peaks by inverse intensity for counting noise.

Only the named element is fitted; all others are assumed natural. Isotopes that
share a nominal mass shift (for example 13C and 2H in one window) cannot be told
apart from a single distribution, and the window must reach the heavy isotope's
mass shift. `ValueError` is raised for an element not in the formula, an
invalid distribution, or a shift beyond the window.

**References:** Kelleher and Masterson, *Am. J. Physiol.* 262, E118 (1992);
Hellerstein and Neese, *Am. J. Physiol.* 263, E988 (1992) (mass isotopomer
distribution analysis).

## Axis length and spacing

`isolen` controls both columns, so the mass axis always matches the intensity
vector length. Override the default spacing through a mass keyword:

```python
distribution = isogen.isodist(
    "AUGC",
    type="RNA",
    isolen=16,
    isotope_spacing=1.0029,
)
```

Formula spacing is editable in the same way:

```python
formula_distribution = isogen.isodist(
    "C6H12O6",
    type="ATOM",
    isotope_spacing=1.0033,
)
```

## Fragment ions and termini

Mass options are forwarded by `isodist`. For example:

```python
b_ion = isogen.isodist("PEPTID", ion_type="b")
triphosphate_rna = isogen.isodist(
    "AUGC",
    type="RNA",
    fiveend="TP",
    threeend="OH",
)
```

For a/b/c ions, pass the N-terminal fragment sequence. For x/y/z ions, pass
the C-terminal fragment sequence. `ion_type` applies only to peptides; named
RNA fragment-ion series are not currently implemented. RNA `threeend` and
`fiveend` values describe supported terminal chemistry for an intact or
manually truncated sequence.

These keyword arguments modify the monoisotopic mass-axis origin. They do not
change the terminal composition used by the FFT, BRAIN, or NN sequence
intensity model.
