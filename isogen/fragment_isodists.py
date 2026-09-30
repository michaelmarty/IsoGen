"""Batch sequence-composition isotope envelopes for peptide fragments."""

from dataclasses import dataclass
from functools import lru_cache
import re

import numpy as np

if __package__:
    from . import mass, protein_mods
    from .isogenwrapper import fft_gen_pep_formula_batch
else:
    import mass
    import protein_mods
    from isogenwrapper import fft_gen_pep_formula_batch


_ELEMENTS = ("C", "H", "N", "O", "S")
if __package__:
    from .protein_composition import RESIDUE_FORMULAS, parse_composition_formula, _tag_composition
else:
    from protein_composition import RESIDUE_FORMULAS, parse_composition_formula, _tag_composition

_RESIDUES = {aa: tuple(parse_composition_formula(RESIDUE_FORMULAS[aa]).get((element, None), 0)
                      for element in _ELEMENTS)
             for aa in protein_mods.CANONICAL_AMINO_ACIDS}
_TERMINI = {
    "a": (-1, 0, 0, -1, 0), "a+1": (-1, 1, 0, -1, 0),
    "b": (0, 0, 0, 0, 0), "c": (0, 3, 1, 0, 0),
    "x": (1, 0, 0, 2, 0), "x+1": (1, 1, 0, 2, 0),
    "y": (0, 2, 0, 1, 0), "y-1": (0, 1, 0, 1, 0),
    "z": (0, -1, -1, 1, 0), "z'": (0, 0, -1, 1, 0),
}
_LABEL = re.compile(r"^([abcxyz])(')?(\d+)([+-]1)?(#\d+)?$")


@dataclass
class FragmentIsotopeBatch:
    labels: tuple
    masses: np.ndarray
    intensities: np.ndarray


def _formula_counts(formula):
    counts = parse_composition_formula(formula)
    if any(element not in _ELEMENTS or isotope is not None
           for element, isotope in counts):
        raise ValueError("fragment modification formula requires natural C, H, N, O, S")
    return np.array([counts.get((element, None), 0) for element in _ELEMENTS], dtype=np.int32)


@lru_cache(maxsize=256)
def _modification_counts(tag, site):
    try:
        counts = _tag_composition(tag, site)
        if any(element not in _ELEMENTS or isotope is not None
               for element, isotope in counts):
            raise ValueError("unsupported fragment isotope element")
        return np.array([counts.get((element, None), 0) for element in _ELEMENTS], dtype=np.int32)
    except ValueError as error:
        raise ValueError("fragment isotope composition is unknown for modification {!r}".format(tag)) from error


def calc_pep_fragment_isodists(sequence, fragmentation_type=None, ion_types=None,
                               *, isolen=128, monoisotopic=True, **kwargs):
    """Return labeled neutral fragment masses and sequence-based FFT intensities.

    Only fragments with an unambiguous C/H/N/O/S composition are supported.
    Mass-only and ambiguously localized modifications raise ``ValueError``.
    """
    if not monoisotopic:
        raise ValueError("fragment isotope envelopes require monoisotopic masses")
    if not isinstance(isolen, int) or not 1 <= isolen <= 128:
        raise ValueError("isolen must be between 1 and 128")
    fragments = mass.calc_pep_fragments(
        sequence, fragmentation_type=fragmentation_type, ion_types=ion_types,
        monoisotopic=True, **kwargs,
    )
    if protein_mods.needs_proforma_parser(sequence):
        plain, modifications, _, _ = protein_mods._parse(sequence)
    else:
        plain, modifications = sequence.upper(), []
    if not plain:
        return FragmentIsotopeBatch((), np.empty(0), np.empty((0, isolen), dtype=np.float32))
    if any(residue not in _RESIDUES for residue in plain):
        raise ValueError("exact fragment isotope composition requires canonical residues")
    residue_counts = np.array([_RESIDUES[residue] for residue in plain], dtype=np.int32)
    prefix = np.vstack((np.zeros(5, dtype=np.int32), np.cumsum(residue_counts, axis=0)))
    total = prefix[-1]
    formulas = []
    lengths = []
    for label in fragments:
        match = _LABEL.fullmatch(label)
        if match is None or match.group(5):
            raise ValueError("ambiguous fragment composition: {!r}".format(label))
        ion_type = match.group(1) + (match.group(2) or "") + (match.group(4) or "")
        length = int(match.group(3))
        if ion_type not in _TERMINI or not 0 < length < len(plain):
            raise ValueError("invalid fragment label: {!r}".format(label))
        is_n_terminal = ion_type[0] in "abc"
        positions = set(range(length) if is_n_terminal else range(len(plain) - length, len(plain)))
        counts = (prefix[length] if is_n_terminal else total - prefix[-length - 1]).copy()
        counts += _TERMINI[ion_type]
        for tag, site, candidates in modifications:
            if tag.strip().startswith("#"):
                continue
            if candidates == "N-term":
                included = is_n_terminal
            elif candidates == "C-term":
                included = not is_n_terminal
            elif candidates is None:
                raise ValueError("unlocalized modification has no exact fragment isotope composition")
            elif candidates <= positions:
                included = True
            elif candidates.isdisjoint(positions):
                included = False
            else:
                raise ValueError("ambiguous fragment modification position")
            if included:
                counts += _modification_counts(tag, site)
        if np.any(counts < 0):
            raise ValueError("fragment has a negative elemental count: {!r}".format(label))
        formulas.append(counts)
        lengths.append(length)
    matrix = fft_gen_pep_formula_batch(np.asarray(formulas, dtype=np.int32).reshape(-1, 5), lengths, isolen)
    return FragmentIsotopeBatch(tuple(fragments), np.asarray(tuple(fragments.values()), dtype=np.float64), matrix)
