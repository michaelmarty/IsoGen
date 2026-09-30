"""Independent isotope reference using Pyteomics enumeration and abundance data.

Enumerate each element separately, then convolve its polynomial. This avoids
the Cartesian explosion of enumerating entire peptide isotopologues. Fixed
isotopes contribute to the mass origin and have no stochastic width.
"""

import numpy as np
from pyteomics import mass


def pyteomics_envelope(natural, fixed=None, isolen=128):
    result = np.zeros(isolen)
    result[0] = 1
    origin = 0.0
    for element, count in natural.items():
        if count < 0:
            raise ValueError("Negative reference atom count")
        if not count:
            continue
        light = min(number for number, (_, abundance) in mass.nist_mass[element].items()
                    if number and abundance > 0)
        origin += count * mass.nist_mass[element][light][0]
        polynomial = np.zeros(isolen)
        for composition, abundance in mass.isotopologues(
                composition={element: count}, report_abundance=True,
                isotope_threshold=1e-15, overall_threshold=0):
            shift = sum((int(name.split("[")[1][:-1]) - light) * number
                        for name, number in composition.items())
            if shift < isolen:
                polynomial[shift] += abundance
        result = np.convolve(result, polynomial)[:isolen]
    for (element, isotope), count in (fixed or {}).items():
        origin += count * mass.nist_mass[element][isotope][0]
    return origin, result / result.max()


def independently_modified(sequence, deltas=(), replacements=(), global_labels=None):
    natural = dict(mass.Composition(sequence=sequence))
    for formula in deltas:
        for element, count in mass.Composition(formula=formula).items():
            natural[element] = natural.get(element, 0) + count
    fixed = {}
    for element, isotope, count in replacements:
        natural[element] -= count
        fixed[(element, isotope)] = fixed.get((element, isotope), 0) + count
    for element, isotope in (global_labels or {}).items():
        fixed[(element, isotope)] = fixed.get((element, isotope), 0) + natural.pop(element, 0)
    return natural, fixed
