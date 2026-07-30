"""IsoGen isotope-distribution calculations for proteins and nucleic acids."""

if __package__:
    from ._version import __version__
    from .isogen import isodist
    from .mass import (
       calc_atom_mass_axis,
       calc_atom_monoisotopic_mass,
       calc_dna_mass,
       calc_dna_mass_axis,
       calc_dna_monoisotopic_mass,
       calc_mass_axis,
       calc_pep_mass,
       calc_pep_mass_axis,
       calc_pep_monoisotopic_mass,
       calc_rna_mass,
       calc_rna_mass_axis,
       calc_rna_monoisotopic_mass,
       gen_mass_axis,
    )
else:
    from _version import __version__
    from isogen import isodist
    from mass import (
       calc_atom_mass_axis,
       calc_atom_monoisotopic_mass,
       calc_dna_mass,
       calc_dna_mass_axis,
       calc_dna_monoisotopic_mass,
       calc_mass_axis,
       calc_pep_mass,
       calc_pep_mass_axis,
       calc_pep_monoisotopic_mass,
       calc_rna_mass,
       calc_rna_mass_axis,
       calc_rna_monoisotopic_mass,
       gen_mass_axis,
    )

__all__ = [
    "__version__",
    "isodist",
    "calc_atom_mass_axis",
    "calc_atom_monoisotopic_mass",
    "calc_dna_mass",
    "calc_dna_mass_axis",
    "calc_dna_monoisotopic_mass",
    "calc_mass_axis",
    "calc_pep_mass",
    "calc_pep_mass_axis",
    "calc_pep_monoisotopic_mass",
    "calc_rna_mass",
    "calc_rna_mass_axis",
    "calc_rna_monoisotopic_mass",
    "gen_mass_axis",
]
