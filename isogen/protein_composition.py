"""Offline ProForma composition resolution; native code receives atoms only."""

from dataclasses import dataclass
from collections import Counter
from functools import lru_cache
import math
import re
import sys

if __package__:
    from . import protein_mods as mods
else:
    import protein_mods as mods


# Residue formulas exclude the terminal water. Keep these independent of the
# legacy native sequence tables so use_modifications=False preserves that path.
RESIDUE_FORMULAS = {
    "A": "C3H5NO", "C": "C3H5NOS", "D": "C4H5NO3", "E": "C5H7NO3",
    "F": "C9H9NO", "G": "C2H3NO", "H": "C6H7N3O", "I": "C6H11NO",
    "K": "C6H12N2O", "L": "C6H11NO", "M": "C5H9NOS", "N": "C4H6N2O2",
    "P": "C5H7NO", "Q": "C5H8N2O2", "R": "C6H12N4O", "S": "C3H5NO2",
    "T": "C4H7NO2", "V": "C5H9NO", "W": "C11H10N2O", "Y": "C9H9NO2",
    "J": "C6H11NO", "U": "C3H5NOSe", "O": "C12H19N3O2", "X": "",
}
_TOKEN = re.compile(
    r"(?:(?:\[(\d+)([A-Z][a-z]?)(-?\d*)\])|(?:(?:\((\d+)\)|(\d+))?([A-Z][a-z]?)))"
    r"(?:\s*\((-?\d+)\)|\s*(-?\d*))"
)
_GLOBAL_ISOTOPE = re.compile(r"(\d+)([A-Z][a-z]?)")
# Molecular subunits in the bundled Unimod formulas. These are residue/delta
# formulas from https://www.unimod.org/xml/unimod.xml (not free sugars).
# Ac means acetate in Unimod, but remains actinium in an elemental Formula:.
_UNIMOD_BRICKS = {
    "Hex": "C6H10O5", "HexNAc": "C8H13NO5", "dHex": "C6H10O4",
    "Pent": "C5H8O4", "HexA": "C6H8O6", "HexN": "C6H11NO4",
    "NeuAc": "C11H17NO8", "NeuGc": "C11H17NO9", "Hep": "C7H12O6",
    "Kdn": "C9H14O8", "Kdo": "C8H12O7", "Ac": "C2H2O", "Me": "CH2",
    "Phos": "HPO3", "Sulf": "SO3", "Water": "H2O",
}


def print_modification_warning(message):
    """Keep warning printouts out of CSV data written to standard output."""
    print("Warning: " + message, file=sys.stderr)


def split_global_isotopes(sequence):
    """Remove isotope replacement headers, preserving fixed modification headers."""
    text = sequence.strip()
    retained = []
    labels = {}
    while text.startswith("<"):
        end = mods._matching(text, 0, "<", ">")
        content = text[1:end]
        if content == "D":
            content = "2H"
        match = _GLOBAL_ISOTOPE.fullmatch(content)
        if match:
            number, element = match.groups()
            if element in labels and labels[element] != int(number):
                raise ValueError("Conflicting global isotope labels for " + element)
            labels[element] = int(number)
        else:
            retained.append(text[:end + 1])
        text = text[end + 1:]
    return "".join(retained) + text, labels


@lru_cache(maxsize=4096)
def _parsed_formula(formula):
    compact = formula.strip()
    if not compact:
        raise ValueError("Formula cannot be empty")
    position = 0
    counts = {}
    for match in _TOKEN.finditer(compact):
        if compact[position:match.start()].strip():
            raise ValueError("Unsupported elemental formula: {!r}".format(formula))
        bracket_number, bracket_element, bracket_count, paren_number, number, element, paren_count, count = match.groups()
        element = bracket_element or element
        isotope = bracket_number or paren_number or number
        isotope = int(isotope) if isotope else None
        if element == "D" and isotope is None:
            element, isotope = "H", 2
        if isotope is not None and not 0 < isotope <= 2147483647:
            raise ValueError("Invalid isotope mass number")
        if element not in mods._ELEMENTS:
            raise ValueError("Unsupported element: " + element)
        count = paren_count if paren_count is not None else count
        if bracket_count:
            if count not in {None, ""}:
                raise ValueError("Isotope count is specified twice")
            count = bracket_count
        count = int(count) if count not in {None, ""} else 1
        key = (element, isotope)
        counts[key] = counts.get(key, 0) + count
        position = match.end()
    if position != len(compact):
        raise ValueError("Unsupported elemental formula: {!r}".format(formula))
    return tuple(sorted(counts.items(), key=lambda item: (item[0][0], item[0][1] or 0)))


def parse_composition_formula(formula):
    """Parse ProForma, Unimod and PSI-MOD signed element/isotope formulas."""
    return dict(_parsed_formula(formula))


def composition_mass(counts, monoisotopic=True):
    if __package__:
        from .mass import atom_masses_monoisotopic
        from .isogenwrapper import isotope_mass
    else:
        from mass import atom_masses_monoisotopic
        from isogenwrapper import isotope_mass

    masses = dict(zip(mods._ELEMENTS, atom_masses_monoisotopic)) if monoisotopic else mods._AVERAGE_ATOM_MASSES
    total = 0.0
    for (element, isotope), count in counts.items():
        if not count:
            continue
        atom_mass = isotope_mass(element, isotope) if isotope else masses.get(element)
        if atom_mass is None:
            raise ValueError("No mass is available for " + element)
        total += atom_mass * count
    return total


def _unimod_composition(formula):
    counts = {}
    for token in formula.split():
        match = re.fullmatch(r"([A-Za-z]+)(?:\((-?\d+)\))?", token)
        if match and match.group(1) in _UNIMOD_BRICKS:
            brick, multiplier = match.groups()
            _add(counts, parse_composition_formula(_UNIMOD_BRICKS[brick]), int(multiplier or 1))
        else:
            _add(counts, parse_composition_formula(token))
    return counts


def _record_composition(record, site):
    if record["cv"] == "RESID":
        choices = record["corrections"]
        if site in choices:
            formula = choices[site]["formula"]
        elif site is None and len(choices) == 1:
            formula = next(iter(choices.values()))["formula"]
        else:
            raise ValueError("Site-dependent RESID composition is unresolved")
    else:
        if not mods._site_matches(record, site):
            raise ValueError("Site-incompatible modification: " + record["accession"])
        formula = record.get("formula")
    if not formula:
        raise ValueError("Modification has no atomic formula")
    # RESID's trailing charge is metadata, not an extra atom or proton.
    formula = re.sub(r"\s+[+-]$", "", formula)
    if record["cv"] == "UNIMOD":
        return _unimod_composition(formula)
    return parse_composition_formula(formula)


@lru_cache(maxsize=4096)
def _lookup_composition(identifier, site, cv=None):
    accessions, names = mods._database()
    accession = identifier.upper().replace("PSI-MOD:", "MOD:")
    if cv:
        records = names[cv].get(identifier.casefold(), [])
    elif accession in accessions:
        records = [accessions[accession]]
    else:
        records = []
        for vocabulary in ("UNIMOD", "MOD", "RESID"):
            records = names[vocabulary].get(identifier.casefold(), [])
            if records:
                break
    compatible = [record for record in records if mods._site_matches(record, site)]
    if not compatible:
        raise ValueError("Unknown or site-incompatible modification: " + identifier)
    formulas = [{key: count for key, count in _record_composition(record, site).items() if count}
                for record in compatible]
    if any(value != formulas[0] for value in formulas[1:]):
        raise ValueError("Ambiguous atomic composition for " + identifier)
    return tuple(formulas[0].items())


def _descriptor_composition(descriptor, site):
    descriptor = mods._AMBIGUITY_SUFFIX.sub("", descriptor.strip()).strip()
    if not descriptor or descriptor.startswith("#"):
        return None
    prefix, separator, value = descriptor.partition(":")
    prefix = prefix.upper()
    if prefix in {"INFO", "OBS"}:
        return None
    if mods._NUMBER.fullmatch(descriptor):
        return None
    if prefix == "FORMULA":
        return parse_composition_formula(value)
    if prefix == "GLYCAN":
        return _glycan_composition(value)
    if prefix in {"GNO", "GNOME", "XLMOD", "XMOD"}:
        raise ValueError(prefix + " composition is not supported")
    if separator and prefix in {"U", "UNIMOD", "M", "MOD", "PSI-MOD", "R", "RESID"}:
        cv = {"U": "UNIMOD", "M": "MOD", "PSI-MOD": "MOD", "R": "RESID"}.get(prefix, prefix)
        if cv == "UNIMOD" and value.isdigit():
            descriptor, cv = "UNIMOD:" + value, None
        elif cv == "MOD" and value.isdigit():
            descriptor, cv = "MOD:" + value.zfill(5), None
        elif cv == "RESID" and value.upper().startswith("AA"):
            descriptor, cv = "RESID:" + value.upper(), None
        else:
            descriptor = value
        return dict(_lookup_composition(descriptor, site, cv))
    try:
        return dict(_lookup_composition(descriptor, site))
    except ValueError as error:
        # Compatibility with the native helper's historical [O1] notation.
        if re.search(r"\d", descriptor):
            try:
                return parse_composition_formula(descriptor)
            except ValueError:
                pass
        raise error


def _glycan_composition(value):
    # Attached monosaccharide residues (free sugar minus H2O).
    formulas = {key: value for key, value in _UNIMOD_BRICKS.items()
                if key not in {"Ac", "Me", "Phos", "Sulf", "Water"}}
    formulas["Fuc"] = "C6H10O4"
    compact = re.sub(r"\s+", "", value)
    pattern = re.compile(r"(" + "|".join(sorted(formulas, key=len, reverse=True)) + r")(\d*)")
    total = {}
    position = 0
    for match in pattern.finditer(compact):
        if match.start() != position:
            raise ValueError("Unsupported glycan composition: " + value)
        name, count = match.groups()
        _add(total, parse_composition_formula(formulas[name]), int(count or 1))
        position = match.end()
    if not compact or position != len(compact):
        raise ValueError("Unsupported glycan composition: " + value)
    return total


def _tag_composition(tag, site):
    if tag.strip().startswith("#"):
        return {}
    values, errors = [], []
    for descriptor in tag.split("|"):
        try:
            value = _descriptor_composition(descriptor, site)
            if value is not None:
                values.append({key: count for key, count in value.items() if count})
        except ValueError as error:
            errors.append(error)
    if not values:
        if errors:
            raise errors[0]
        if all(part.strip().lower().startswith(("info:", "obs:", "#")) for part in tag.split("|")):
            return {}
        raise ValueError("A mass shift alone does not specify atomic composition")
    if any(value != values[0] for value in values[1:]):
        raise ValueError("Conflicting atomic compositions in tag")
    mono = composition_mass(values[0])
    for descriptor in tag.split("|"):
        cleaned = mods._AMBIGUITY_SUFFIX.sub("", descriptor).strip()
        if mods._NUMBER.fullmatch(cleaned) and not math.isclose(mono, float(cleaned), abs_tol=0.01):
            raise ValueError("Conflicting masses in modification tag: " + tag)
    return values[0]


def _add(target, delta, multiplier=1):
    for key, count in delta.items():
        target[key] = target.get(key, 0) + multiplier * count


def _terminal_delta(ion_type, labels=None):
    """Terminal composition correction relative to the intact H2O formula."""
    if __package__:
        from .mass import _normalize_pep_ion_type
    else:
        from mass import _normalize_pep_ion_type
    formulas = {"h2o": "H2O", "a": "C-1O-1", "a+1": "C-1O-1H",
                "b": "", "c": "NH3", "x": "CO2", "x+1": "CO2H",
                "y": "H2O", "y-1": "HO", "z": "H-1N-1O", "z'": "N-1O"}
    formula = formulas[_normalize_pep_ion_type(ion_type)]
    delta = parse_composition_formula(formula) if formula else {}
    _add(delta, parse_composition_formula("H2O"), -1)
    if labels:
        return {(element, labels.get(element)): count for (element, _), count in delta.items()}
    return delta


@dataclass
class ProteinComposition:
    sequence: str
    counts: dict
    mass: float
    changed: bool


def resolve_proforma_composition(sequence, use_modifications=True, ion_type="H2O", warn=True):
    """Resolve intact composition, skipping unsupported annotations individually.

    Mass-only annotations still contribute their known mass to the axis. Unknown
    descriptors contribute neither mass nor atoms; malformed syntax still raises.
    """
    if __package__:
        from . import mass as base_mass
    else:
        import mass as base_mass
    text, global_labels = split_global_isotopes(sequence)
    plain, annotations, _, terminal = mods._parse(text)
    if terminal and str(ion_type).upper() != "H2O":
        raise ValueError("Explicit terminal modifications cannot be combined with fragment ion_type")
    counts = parse_composition_formula("H2O")
    residue_masses = dict(base_mass.aa_masses_monoisotopic)
    residue_masses.update(J=residue_masses["I"], U=150.953634, O=237.147727,
                          B=(residue_masses["D"] + residue_masses["N"]) / 2,
                          Z=(residue_masses["E"] + residue_masses["Q"]) / 2, X=0.0)
    # General-element envelopes start at the lightest isotope, whereas legacy
    # peptide/CV masses may use the most abundant isotope (notably selenium).
    residue_masses["U"] = composition_mass(parse_composition_formula(RESIDUE_FORMULAS["U"]))
    residue_counts = Counter(plain)
    total_mass = sum(residue_masses[aa] * count for aa, count in residue_counts.items())
    total_mass += base_mass.get_pep_ion_mass_shift(ion_type, monoisotopic=True)
    for aa, count in residue_counts.items():
        if aa in {"B", "Z"}:
            if use_modifications and warn:
                print_modification_warning("Ignoring ambiguous residue {} for isotope prediction".format(aa))
            continue
        formula = RESIDUE_FORMULAS[aa]
        if formula:
            _add(counts, parse_composition_formula(formula), count)
    terminal_delta = _terminal_delta(ion_type)
    _add(counts, terminal_delta)
    if any(count < 0 for count in counts.values()):
        raise ValueError("Terminal chemistry produces negative atom counts")
    changed = any(aa in "JOU" for aa in plain) or str(ion_type).lower() != "h2o"
    for tag, site, _ in annotations:
        delta = None
        problem = None
        try:
            delta = _tag_composition(tag, site)
        except ValueError as error:
            problem = str(error)
        try:
            tag_mass = mods._tag_mass(tag, site, True)
        except ValueError:
            tag_mass = composition_mass(delta) if delta is not None else 0.0
        if delta is not None:
            tag_mass = composition_mass(delta)
        total_mass += tag_mass
        if delta is not None:
            candidate = counts.copy()
            _add(candidate, delta)
            if any(count < 0 for count in candidate.values()):
                problem = "Modification would produce negative atom counts"
            elif any(count > 2147483647 for count in candidate.values()):
                problem = "Modification atom count exceeds the native integer range"
            else:
                counts = candidate
                changed |= any(delta.values())
        if use_modifications and problem and warn:
            print_modification_warning("Ignoring modification {!r} for isotope prediction: {}".format(tag, problem))
    for element, isotope in global_labels.items():
        try:
            labeled_mass = composition_mass({(element, isotope): 1})
            natural_mass = composition_mass({(element, None): 1})
        except ValueError as error:
            if use_modifications and warn:
                print_modification_warning("Ignoring global isotope label: " + str(error))
            continue
        count = counts.get((element, None), 0)
        total_mass += count * (labeled_mass - natural_mass)
        counts[(element, None)] = 0
        counts[(element, isotope)] = counts.get((element, isotope), 0) + count
        changed |= count > 0
    return ProteinComposition(plain, {key: count for key, count in counts.items() if count}, total_mass, changed)
