import numpy as np
import re
import os

aa_masses = {'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766,
             'G': 57.0519, 'H': 137.1411, 'I': 113.1594, 'K': 128.1741, 'L': 113.1594,
             'M': 131.1926, 'N': 114.1038, 'P': 97.1167, 'Q': 128.1307, 'R': 156.1875,
             'S': 87.0782, 'T': 101.1051, 'V': 99.1326, 'W': 186.2132, 'Y': 163.1760}

aa_masses_monoisotopic = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,
                          'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,
                          'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
                          'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333}

rna_masses = {'A': 329.2, 'U': 306.2, 'C': 305.2, 'G': 345.2, 'T': 306.2}

dna_masses = {'A': 313.2, 'T': 304.2, 'C': 289.2, 'G': 329.2, 'U': 304.2, }

rna_masses_monoisotopic = {'A': 329.05252, 'U': 306.02530, 'C': 305.04129, 'G': 345.04744,
                           'T': 306.02530}

dna_masses_monoisotopic = {'A': 313.05761, 'T': 304.04604, 'C': 289.04637, 'G': 329.05252,
                           'U': 304.04604}

mass_water = 18.0153
mass_OH = 17.008
mass_O = 15.9994
mass_HPO4 = 95.9793
mass_H = 1.00794
mass_proton = 1.00727647
mass_CO = 28.0101
mass_NH3 = 17.03052
mass_CO2 = 44.0095

# Monoisotopic terminal-group masses.
mass_water_monoisotopic = 18.010564684
mass_OH_monoisotopic = 17.002739652
mass_O_monoisotopic = 15.994914620
mass_HPO4_monoisotopic = 95.961245509
mass_H_monoisotopic = 1.007825032
mass_CO_monoisotopic = 27.994914620
mass_NH3_monoisotopic = 17.026549101
mass_CO2_monoisotopic = 43.989829239

# Neutral terminal-group shifts relative to the sum of amino-acid residue
# masses. a/b/c are N-terminal fragments; x/y/z are C-terminal fragments.
pep_ion_mass_shifts = {
    "H2O": mass_water,
    "A": -mass_CO,
    "B": 0.0,
    "C": mass_NH3,
    "X": mass_CO2,
    "Y": mass_water,
    "Z": mass_water - mass_NH3,
}

pep_ion_mass_shifts_monoisotopic = {
    "H2O": mass_water_monoisotopic,
    "A": -mass_CO_monoisotopic,
    "B": 0.0,
    "C": mass_NH3_monoisotopic,
    "X": mass_CO2_monoisotopic,
    "Y": mass_water_monoisotopic,
    "Z": mass_water_monoisotopic - mass_NH3_monoisotopic,
}


def get_aa_mass(letter):
    """Return the average residue mass for a one-letter amino-acid code.

    Whitespace and unknown codes contribute zero mass; unknown codes are also
    reported to standard output.
    """
    if letter == " " or letter == "\t" or letter == "\n":
        return 0

    try:
        return aa_masses[letter]
    except Exception as exception:
        print("Bad Amino Acid Code:", letter)
        return 0


def get_rna_mass(letter):
    """Return the average RNA residue mass for a nucleotide code.

    ``T`` is accepted as uracil. Unknown codes contribute zero mass.
    """
    if letter == "T":
        print("Assuming T means U")

    try:
        return rna_masses[letter]
    except Exception as exception:
        print("Bad RNA Code:", letter)
        return 0


def get_dna_mass(letter):
    """Return the average DNA residue mass for a nucleotide code.

    ``U`` is accepted as thymine. Unknown codes contribute zero mass.
    """
    try:
        return dna_masses[letter]
    except Exception as exception:
        print("Bad DNA Code:", letter)
        return 0


def _get_monoisotopic_mass(letter, masses, molecule_name):
    """Look up a monoisotopic residue mass with shared error handling."""
    if letter == " " or letter == "\t" or letter == "\n":
        return 0

    try:
        return masses[letter]
    except KeyError:
        print("Bad " + molecule_name + " Code:", letter)
        return 0


def get_aa_monoisotopic_mass(letter):
    """Return the monoisotopic residue mass for an amino-acid code."""
    return _get_monoisotopic_mass(letter, aa_masses_monoisotopic, "Amino Acid")


def get_rna_monoisotopic_mass(letter):
    """Return an RNA residue's monoisotopic mass, treating ``T`` as ``U``."""
    if letter == "T":
        print("Assuming T means U")
    return _get_monoisotopic_mass(letter, rna_masses_monoisotopic, "RNA")


def get_dna_monoisotopic_mass(letter):
    """Return a DNA residue's monoisotopic mass, treating ``U`` as ``T``."""
    return _get_monoisotopic_mass(letter, dna_masses_monoisotopic, "DNA")


def get_pep_ion_mass_shift(ion_type="H2O", monoisotopic=False):
    """Return the neutral terminal-group mass shift for a protein fragment.

    ``H2O`` represents an intact protein. The supplied sequence should be the
    N-terminal fragment for a/b/c ions or the C-terminal fragment for x/y/z
    ions.

    Args:
        ion_type: ``H2O`` for an intact protein, or ``a``, ``b``, ``c``,
            ``x``, ``y``, or ``z`` for a fragment ion.
        monoisotopic: Use monoisotopic rather than average atomic masses.

    Returns:
        Signed terminal-group mass shift in daltons.

    Raises:
        TypeError: If ``ion_type`` is not a string.
        ValueError: If the ion type is unsupported.
    """
    if not isinstance(ion_type, str):
        raise TypeError("ion_type must be a string")

    normalized_type = ion_type.upper()
    shifts = pep_ion_mass_shifts_monoisotopic if monoisotopic else pep_ion_mass_shifts
    try:
        return shifts[normalized_type]
    except KeyError as exception:
        choices = ", ".join(pep_ion_mass_shifts)
        raise ValueError("Unknown ion_type {!r}; expected one of: {}".format(ion_type, choices)) from exception


def calc_pep_mass(sequence, allow_float=True, remove_nan=True, all_cyst_ox=False, pyroglu=False, round_to=2,
                  ion_type="H2O"):
    """Calculate an average protein or protein-fragment mass.

    Args:
        sequence: Amino-acid sequence or an existing numeric mass.
        allow_float: Interpret numeric input as an already calculated mass.
        remove_nan: Return zero for the string ``"nan"``.
        all_cyst_ox: Remove one hydrogen mass per cysteine.
        pyroglu: Apply an N-terminal pyroglutamate loss for glutamate or
            glutamine.
        round_to: Number of decimal places in the returned value.
        ion_type: ``H2O`` for an intact protein, or a/b/c/x/y/z for a supplied
            N- or C-terminal fragment sequence.

    Returns:
        Average neutral mass in daltons.
    """
    is_sequence = isinstance(sequence, str)
    if all_cyst_ox and is_sequence:
        # Count number of c in sequence
        c = sequence.lower().count("c")
        # Multiply by -1 * mass of H
        modmass = c * (-1 * mass_H)
    else:
        modmass = 0

    if remove_nan and is_sequence and sequence.lower() == "nan":
        return 0.0

    if allow_float:
        try:
            mass = float(sequence)
        except (TypeError, ValueError):
            seq = sequence.upper()
            mass = np.sum([get_aa_mass(s) for s in seq])
            mass += get_pep_ion_mass_shift(ion_type)
    else:
        seq = sequence.upper()
        mass = np.sum([get_aa_mass(s) for s in seq])
        mass += get_pep_ion_mass_shift(ion_type)
    # print(sequence, mass)
    # Look for pyroglutamate mod if set
    if pyroglu and is_sequence and len(sequence) > 0:
        if sequence[0].upper() == "E":
            modmass -= mass_water
        if sequence[0].upper() == "Q":
            modmass -= mass_OH

    massoutput = np.round(mass + modmass, round_to)
    return massoutput


def calc_pep_monoisotopic_mass(sequence, allow_float=True, remove_nan=True, all_cyst_ox=False, pyroglu=False,
                               ion_type="H2O"):
    """Calculate a monoisotopic protein or protein-fragment mass.

    Args:
        sequence: Amino-acid sequence or an existing numeric mass.
        allow_float: Interpret numeric input as an already calculated mass.
        remove_nan: Return zero for the string ``"nan"``.
        all_cyst_ox: Remove one monoisotopic hydrogen mass per cysteine.
        pyroglu: Apply an N-terminal pyroglutamate loss for glutamate or
            glutamine.
        ion_type: ``H2O`` for an intact protein, or a/b/c/x/y/z for a supplied
            N- or C-terminal fragment sequence.

    Returns:
        Monoisotopic neutral mass in daltons.
    """
    is_sequence = isinstance(sequence, str)
    if remove_nan and is_sequence and sequence.lower() == "nan":
        return 0.0

    if allow_float:
        try:
            mass = float(sequence)
        except (TypeError, ValueError):
            seq = sequence.upper()
            mass = np.sum([get_aa_monoisotopic_mass(s) for s in seq])
            mass += get_pep_ion_mass_shift(ion_type, monoisotopic=True)
    else:
        seq = sequence.upper()
        mass = np.sum([get_aa_monoisotopic_mass(s) for s in seq])
        mass += get_pep_ion_mass_shift(ion_type, monoisotopic=True)

    modmass = 0.0
    if all_cyst_ox and is_sequence:
        modmass -= sequence.lower().count("c") * mass_H_monoisotopic

    if pyroglu and is_sequence and len(sequence) > 0:
        if sequence[0].upper() == "E":
            modmass -= mass_water_monoisotopic
        elif sequence[0].upper() == "Q":
            modmass -= mass_OH_monoisotopic

    return float(mass + modmass)


def calc_rna_mass(sequence, threeend="OH", fiveend="MP"):
    """Calculate the average neutral mass of an RNA sequence.

    Args:
        sequence: RNA sequence; ``T`` is treated as ``U``.
        threeend: Three-prime terminus, currently ``OH`` or no adjustment.
        fiveend: Five-prime terminus: ``OH``, monophosphate (``MP``), or
            triphosphate (``TP``).

    Returns:
        Average neutral mass in daltons.
    """
    seq = sequence.upper()
    mass = np.sum([get_rna_mass(s) for s in seq])
    if threeend == "OH":
        mass += mass_OH

    if fiveend == "OH":
        mass -= mass_HPO4
        mass += mass_OH
    elif fiveend == "MP":
        mass += mass_H
    elif fiveend == "TP":
        mass += mass_HPO4 + mass_HPO4 - mass_O - mass_O + mass_H

    return float(mass)


def calc_rna_monoisotopic_mass(sequence, threeend="OH", fiveend="MP"):
    """Calculate the monoisotopic neutral mass of an RNA sequence.

    Args:
        sequence: RNA sequence; ``T`` is treated as ``U``.
        threeend: Three-prime terminus, currently ``OH`` or no adjustment.
        fiveend: Five-prime terminus: ``OH``, monophosphate (``MP``), or
            triphosphate (``TP``).

    Returns:
        Monoisotopic neutral mass in daltons.
    """
    seq = sequence.upper()
    mass = np.sum([get_rna_monoisotopic_mass(s) for s in seq])
    if threeend == "OH":
        mass += mass_OH_monoisotopic

    if fiveend == "OH":
        mass -= mass_HPO4_monoisotopic
        mass += mass_OH_monoisotopic
    elif fiveend == "MP":
        mass += mass_H_monoisotopic
    elif fiveend == "TP":
        mass += (2 * mass_HPO4_monoisotopic) - (2 * mass_O_monoisotopic) + mass_H_monoisotopic

    return float(mass)


def calc_dna_mass(sequence, threeend="OH", fiveend="MP"):
    """Calculate the average neutral mass of a DNA sequence.

    Args:
        sequence: DNA sequence; ``U`` is treated as ``T``.
        threeend: Three-prime terminus, currently ``OH`` or no adjustment.
        fiveend: Five-prime terminus: ``OH``, monophosphate (``MP``), or
            triphosphate (``TP``).

    Returns:
        Average neutral mass in daltons.
    """
    seq = sequence.upper()
    mass = np.sum([get_dna_mass(s) for s in seq])
    if threeend == "OH":
        mass += mass_OH

    if fiveend == "OH":
        mass -= mass_HPO4
        mass += mass_OH
    elif fiveend == "MP":
        mass += mass_H
    elif fiveend == "TP":
        mass += mass_HPO4 + mass_HPO4 - mass_O - mass_O + mass_H

    return float(mass)


def calc_dna_monoisotopic_mass(sequence, threeend="OH", fiveend="MP"):
    """Calculate the monoisotopic neutral mass of a DNA sequence.

    Args:
        sequence: DNA sequence; ``U`` is treated as ``T``.
        threeend: Three-prime terminus, currently ``OH`` or no adjustment.
        fiveend: Five-prime terminus: ``OH``, monophosphate (``MP``), or
            triphosphate (``TP``).

    Returns:
        Monoisotopic neutral mass in daltons.
    """
    seq = sequence.upper()
    mass = np.sum([get_dna_monoisotopic_mass(s) for s in seq])
    if threeend == "OH":
        mass += mass_OH_monoisotopic

    if fiveend == "OH":
        mass -= mass_HPO4_monoisotopic
        mass += mass_OH_monoisotopic
    elif fiveend == "MP":
        mass += mass_H_monoisotopic
    elif fiveend == "TP":
        mass += (2 * mass_HPO4_monoisotopic) - (2 * mass_O_monoisotopic) + mass_H_monoisotopic

    return float(mass)


def calc_mass_axis(monoisotopic_mass, isolen=128, isotope_spacing=1.0033):
    """Create an evenly spaced mass axis.

    Args:
        monoisotopic_mass: First mass-axis value.
        isolen: Number of values to generate.
        isotope_spacing: Mass difference between adjacent isotope positions.

    Returns:
        A one-dimensional float NumPy array.

    Raises:
        TypeError: If ``isolen`` is not an integer.
        ValueError: If ``isolen`` is negative.
    """
    if not isinstance(isolen, (int, np.integer)):
        raise TypeError("isolen must be an integer")
    if isolen < 0:
        raise ValueError("isolen cannot be negative")

    return float(monoisotopic_mass) + np.arange(isolen, dtype=float) * float(isotope_spacing)


def calc_pep_mass_axis(sequence, isolen=128, isotope_spacing=1.0033, ion_type="H2O", **mass_kwargs):
    """Create a protein or fragment mass axis from its sequence.

    Args:
        sequence: Protein or fragment sequence.
        isolen: Number of values to generate.
        isotope_spacing: Difference between adjacent isotope masses.
        ion_type: Intact or fragment-ion terminal composition.
        **mass_kwargs: Additional options for
            :func:`calc_pep_monoisotopic_mass`.

    Returns:
        A one-dimensional float NumPy array.
    """
    monoisotopic_mass = calc_pep_monoisotopic_mass(sequence, ion_type=ion_type, **mass_kwargs)
    return calc_mass_axis(monoisotopic_mass, isolen, isotope_spacing=isotope_spacing)


def calc_rna_mass_axis(sequence, isolen=128, isotope_spacing=1.0027, threeend="OH", fiveend="MP"):
    """Create an RNA mass axis from a sequence and terminal chemistry.

    Args:
        sequence: RNA sequence.
        isolen: Number of values to generate.
        isotope_spacing: Difference between adjacent isotope masses.
        threeend: Three-prime terminal chemistry.
        fiveend: Five-prime terminal chemistry.

    Returns:
        A one-dimensional float NumPy array beginning at the RNA
        monoisotopic mass.
    """
    monoisotopic_mass = calc_rna_monoisotopic_mass(sequence, threeend=threeend, fiveend=fiveend)
    return calc_mass_axis(monoisotopic_mass, isolen, isotope_spacing=isotope_spacing)


def calc_dna_mass_axis(sequence, isolen=128, isotope_spacing=1.0027, threeend="OH", fiveend="MP"):
    """Create a DNA mass axis from a sequence and terminal chemistry.

    Args:
        sequence: DNA sequence.
        isolen: Number of values to generate.
        isotope_spacing: Difference between adjacent isotope masses.
        threeend: Three-prime terminal chemistry.
        fiveend: Five-prime terminal chemistry.

    Returns:
        A one-dimensional float NumPy array beginning at the DNA
        monoisotopic mass.
    """
    monoisotopic_mass = calc_dna_monoisotopic_mass(sequence, threeend=threeend, fiveend=fiveend)
    return calc_mass_axis(monoisotopic_mass, isolen, isotope_spacing=isotope_spacing)


def gen_mass_axis(input, type="PEPTIDE", isolen=128, isotope_spacing=None, **mass_kwargs):
    """Create a mass axis from a numeric mass or biopolymer sequence.

    Args:
        input: Numeric monoisotopic mass or sequence string.
        type: ``PEPTIDE``, ``RNA``, or ``DNA``.
        isolen: Number of values to generate.
        isotope_spacing: Optional spacing override. Defaults to 1.0033 Da for
            proteins and 1.0027 Da for nucleic acids.
        **mass_kwargs: Sequence mass options such as ``ion_type``,
            ``threeend``, or ``fiveend``.

    Returns:
        A one-dimensional float NumPy array.
    """
    if not isinstance(input, str):
        if isotope_spacing is None:
            isotope_spacing = 1.0027 if type in ("RNA", "DNA") else 1.0033
        return calc_mass_axis(input, isolen, isotope_spacing=isotope_spacing)

    if type == "PEPTIDE":
        if isotope_spacing is None:
            isotope_spacing = 1.0033
        return calc_pep_mass_axis(input, isolen, isotope_spacing=isotope_spacing, **mass_kwargs)
    elif type == "RNA":
        if isotope_spacing is None:
            isotope_spacing = 1.0027
        return calc_rna_mass_axis(input, isolen, isotope_spacing=isotope_spacing, **mass_kwargs)
    elif type == "DNA":
        if isotope_spacing is None:
            isotope_spacing = 1.0027
        return calc_dna_mass_axis(input, isolen, isotope_spacing=isotope_spacing, **mass_kwargs)
    else:
        raise ValueError("Unknown type for mass axis generation: {}".format(type))


def read_fasta(path):
    """Read a FASTA file into a mapping of identifiers to sequences.

    Args:
        path: Path to a FASTA-formatted text file.

    Returns:
        A dictionary mapping the first token of each header to its concatenated
        sequence.
    """
    f = open(path, 'r')
    lines = f.readlines()

    hre = re.compile(r'>(\S+)')
    lre = re.compile(r'^(\S+)$')

    gene = {}

    for line in lines:
        outh = hre.search(line)
        if outh:
            id = outh.group(1)
        else:
            outl = lre.search(line)
            if id in gene.keys():
                gene[id] += outl.group(1)
            else:
                gene[id] = outl.group(1)
    return gene


if __name__ == "__main__":
    #print(mass_HPO4 + mass_HPO4 - mass_O - mass_O + mass_H + mass_H)
    os.chdir("..\\..\\Scripts\\old\\Jessica")
    file = "ecoli.fasta"
    genes = read_fasta(file)
    print(genes)
    exit()
    oligo = "aacauucaACgcugucggugAgu"
    mass = calc_rna_mass(oligo)
    print(mass)
