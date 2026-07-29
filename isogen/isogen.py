import numpy as np

if __package__:
    from . import isogenwrapper as wrapper
    from . import mass
else:
    import isogenwrapper as wrapper
    import mass


def isodist(input, type="PEPTIDE", isolen=128, method="FFT", **mass_kwargs):
    """
    Generate an isotope distribution for the given input.

    Args:
        input: The input data for which to generate the isotope distribution.
        type: The input type: ``PEPTIDE``, ``RNA``, or ``DNA``. DNA uses the
            RNA intensity model with a DNA-specific mass axis.
        isolen: The length of the isotope distribution (default is 128).
        method: The method to use for generating the isotope distribution (default is "FFT").
        **mass_kwargs: Additional options passed to mass-axis generation, such
            as ``ion_type`` for peptide fragment ions.

    Returns:
        The generated isotope distribution.
    """
    int_dist = wrapper.gen_isodist(input, type=type, isolen=isolen, method=method)
    mass_axis = mass.gen_mass_axis(input, type=type, isolen=isolen, **mass_kwargs)
    output = np.transpose(np.vstack((mass_axis, int_dist)))
    return output

if __name__ == "__main__":
    print("Launching IsoGen")
    test_mass = 10000
    test_pep_seq = "ACDEFGHIKLMNPQRSTVWY"
    test_rna_seq = "AUGCAGUACGUA"

    print(isodist(10000, type="PEPTIDE", isolen=128, method="FFT"))
