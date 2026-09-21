"""Mass calculation for modified proteins written in ProForma notation."""

from functools import lru_cache
import gzip
import json
import math
from pathlib import Path
import re


CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
PROFORMA_AMINO_ACIDS = CANONICAL_AMINO_ACIDS | frozenset("BJOUXZ")

_ELEMENTS = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
    "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti",
    "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs",
    "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
    "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
    "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
)

# Standard atomic weights needed by ordinary biomolecular formulas. Exact
# isotope-zero masses are taken from mass.atom_masses_monoisotopic instead.
_AVERAGE_ATOM_MASSES = {
    "H": 1.00794, "B": 10.811, "C": 12.0107, "N": 14.0067,
    "O": 15.9994, "F": 18.9984032, "Na": 22.9897693, "Mg": 24.305,
    "Al": 26.9815385, "Si": 28.0855, "P": 30.973762, "S": 32.065,
    "Cl": 35.453, "K": 39.0983, "Ca": 40.078, "Fe": 55.845,
    "Se": 78.96, "Br": 79.904, "I": 126.90447,
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(?:\((-?\d+)\)|(-?\d*))")
_NUMBER = re.compile(r"^[+-](?:\d+(?:\.\d*)?|\.\d+)$")
_AMBIGUITY_SUFFIX = re.compile(r"#([A-Za-z0-9_.-]+)(?:\([^)]*\))?")


def needs_proforma_parser(sequence):
    """Return whether a protein string needs the modified-sequence path."""
    if not isinstance(sequence, str):
        return False
    return any(character in "[]{}<>()?-/+|" for character in sequence) or any(
        character.isalpha() and character.upper() not in CANONICAL_AMINO_ACIDS
        for character in sequence
    )


def _matching(text, start, opening, closing):
    depth = 0
    for index in range(start, len(text)):
        if text[index] == opening:
            depth += 1
        elif text[index] == closing:
            depth -= 1
            if depth == 0:
                return index
    raise ValueError("Unclosed {!r} in ProForma sequence".format(opening))


def _outside_character(text, target):
    stack = []
    pairs = {")": "(", "]": "[", "}": "{", ">": "<"}
    for character in text:
        if character in "([{<":
            stack.append(character)
        elif character in pairs:
            if not stack or stack.pop() != pairs[character]:
                raise ValueError("Unbalanced ProForma delimiters")
        elif character == target and not stack:
            return True
    if stack:
        raise ValueError("Unbalanced ProForma delimiters")
    return False


def _parse_core(text):
    sequence = []
    modifications = []
    x_gaps = []
    has_terminal_modification = False
    index = 0
    while index < len(text):
        character = text[index]
        if character.isalpha() and len(character) == 1:
            residue = character.upper()
            if residue not in PROFORMA_AMINO_ACIDS:
                raise ValueError("Unknown amino-acid code: {!r}".format(character))
            sequence.append(residue)
            index += 1
            tags = []
            while index < len(text) and text[index] == "[":
                end = _matching(text, index, "[", "]")
                tag = text[index + 1:end]
                modifications.append((tag, residue))
                tags.append(tag)
                index = end + 1
            if residue == "X":
                x_gaps.append(tags)
            continue
        if character == "(":
            end = _matching(text, index, "(", ")")
            region_sequence, region_mods, region_x, region_terminal = _parse_core(
                text[index + 1:end]
            )
            sequence.extend(region_sequence)
            modifications.extend(region_mods)
            x_gaps.extend(region_x)
            has_terminal_modification |= region_terminal
            index = end + 1
            while index < len(text) and text[index] == "[":
                tag_end = _matching(text, index, "[", "]")
                modifications.append((text[index + 1:tag_end], None))
                index = tag_end + 1
            continue
        if character == "-" and index + 1 < len(text) and text[index + 1] == "[":
            has_terminal_modification = True
            index += 1
            while index < len(text) and text[index] == "[":
                end = _matching(text, index, "[", "]")
                modifications.append((text[index + 1:end], "C-term"))
                index = end + 1
            if index != len(text):
                raise ValueError("Unexpected text after a C-terminal modification")
            continue
        raise ValueError(
            "Unsupported or malformed ProForma syntax near {!r}".format(text[index:])
        )
    return sequence, modifications, x_gaps, has_terminal_modification


def _parse(sequence):
    if not isinstance(sequence, str):
        raise TypeError("ProForma sequence must be a string")
    text = sequence.strip()
    if not text:
        raise ValueError("ProForma sequence cannot be empty")
    if _outside_character(text, "+"):
        raise ValueError("Chimeric and multi-peptidoform expressions are not supported")
    if "//" in text or "#XL" in text.upper() or "#BRANCH" in text.upper():
        raise ValueError("Cross-linked and branched peptides are not supported")

    charge = re.search(r"/[-+]?\d+$", text)
    if charge:
        text = text[:charge.start()]

    modifications = []
    fixed = []
    terminal = False

    while text.startswith("<"):
        end = _matching(text, 0, "<", ">")
        content = text[1:end]
        match = re.fullmatch(r"\[(.*)\]@(.+)", content)
        if not match:
            raise ValueError(
                "Global isotope replacement is not yet supported; expected <[mod]@site>"
            )
        fixed.append((match.group(1), match.group(2)))
        text = text[end + 1:]

    while text.startswith("{"):
        end = _matching(text, 0, "{", "}")
        modifications.append((text[1:end], None))
        text = text[end + 1:]

    while text.startswith("["):
        end = _matching(text, 0, "[", "]")
        remainder = text[end + 1:]
        occurrence = re.match(r"\^(\d+)\?", remainder)
        if remainder.startswith("?") or occurrence:
            count = int(occurrence.group(1)) if occurrence else 1
            if count < 1:
                raise ValueError("Unlocalized modification count must be positive")
            modifications.extend((text[1:end], None) for _ in range(count))
            text = remainder[occurrence.end():] if occurrence else remainder[1:]
        elif remainder.startswith("-"):
            modifications.append((text[1:end], "N-term"))
            terminal = True
            text = text[end + 2:]
        else:
            break

    core_sequence, core_mods, x_gaps, core_terminal = _parse_core(text)
    modifications.extend(core_mods)
    terminal |= core_terminal
    base_sequence = "".join(core_sequence)
    if not base_sequence:
        raise ValueError("ProForma input does not contain an amino-acid sequence")

    for tag, selector_text in fixed:
        selectors = [part.strip() for part in selector_text.split(",")]
        for selector in selectors:
            normalized = selector.lower().replace("protein ", "")
            if len(selector) == 1 and selector.upper() in PROFORMA_AMINO_ACIDS:
                for _ in range(base_sequence.count(selector.upper())):
                    modifications.append((tag, selector.upper()))
            elif normalized in {"n-term", "c-term"}:
                modifications.append((tag, normalized[0].upper() + normalized[1:]))
                terminal = True
            else:
                raise ValueError("Unsupported global modification selector: {!r}".format(selector))

    return base_sequence, modifications, x_gaps, terminal


def strip_proforma(sequence):
    """Return the unannotated amino-acid sequence from a ProForma string."""
    if not needs_proforma_parser(sequence):
        return sequence.upper()
    return _parse(sequence)[0]


@lru_cache(maxsize=1)
def _database():
    path = Path(__file__).with_name("resources") / "protein_modifications.json.gz"
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        records = json.load(handle)["records"]
    accessions = {record["accession"].upper(): record for record in records}
    names = {"UNIMOD": {}, "MOD": {}, "RESID": {}}
    for record in records:
        index = names[record["cv"]]
        for name in [record["name"], *record.get("synonyms", ())]:
            index.setdefault(name.casefold(), []).append(record)
    return accessions, names


def _site_matches(record, site):
    if site is None or record["cv"] == "RESID":
        return True
    sites = record.get("sites") or []
    if not sites or "X" in sites:
        return True
    if site in {"N-term", "C-term"}:
        return any(value.lower().endswith(site.lower()) for value in sites)
    return site in sites


def _record_mass(record, site, monoisotopic):
    key = "mono" if monoisotopic else "average"
    if record["cv"] == "RESID":
        corrections = record["corrections"]
        if site in corrections:
            return corrections[site][key]
        if site is None and len(corrections) == 1:
            return next(iter(corrections.values()))[key]
        choices = ", ".join(sorted(corrections))
        raise ValueError(
            "{} is not defined for site {!r}; expected one of {}".format(
                record["accession"], site, choices
            )
        )
    if not _site_matches(record, site):
        raise ValueError(
            "{} is not valid at site {!r}".format(record["accession"], site)
        )
    value = record.get(key)
    if value is None:
        value = record.get("mono")
    return float(value)


def _lookup(identifier, site, monoisotopic, cv=None):
    accessions, names = _database()
    if cv:
        cv = cv.upper()
        records = names[cv].get(identifier.casefold(), [])
    else:
        accession = identifier.upper().replace("PSI-MOD:", "MOD:")
        if accession in accessions:
            return _record_mass(accessions[accession], site, monoisotopic)
        records = []
        for vocabulary in ("UNIMOD", "MOD", "RESID"):
            records.extend(names[vocabulary].get(identifier.casefold(), []))
            if records:
                break
    compatible = [record for record in records if _site_matches(record, site)]
    if not compatible:
        raise ValueError("Unknown or site-incompatible modification: {!r}".format(identifier))
    masses = [_record_mass(record, site, monoisotopic) for record in compatible]
    if not all(math.isclose(masses[0], value, abs_tol=1e-6) for value in masses[1:]):
        raise ValueError("Ambiguous modification name: {!r}".format(identifier))
    return masses[0]


def _formula_mass(formula, monoisotopic):
    compact = re.sub(r"\s+", "", formula)
    position = 0
    total = 0.0
    if not compact:
        raise ValueError("Formula cannot be empty")
    if monoisotopic:
        if __package__:
            from .mass import atom_masses_monoisotopic
        else:
            from mass import atom_masses_monoisotopic
        masses = dict(zip(_ELEMENTS, atom_masses_monoisotopic))
    else:
        masses = _AVERAGE_ATOM_MASSES
    for match in _FORMULA_TOKEN.finditer(compact):
        if match.start() != position:
            raise ValueError("Unsupported elemental formula: {!r}".format(formula))
        element = match.group(1)
        count_text = match.group(2) if match.group(2) is not None else match.group(3)
        count = int(count_text) if count_text not in {None, ""} else 1
        try:
            total += masses[element] * count
        except KeyError as exception:
            kind = "monoisotopic" if monoisotopic else "average"
            raise ValueError("No {} mass is available for {}".format(kind, element)) from exception
        position = match.end()
    if position != len(compact):
        raise ValueError("Unsupported elemental formula: {!r}".format(formula))
    return total


def _descriptor_mass(descriptor, site, monoisotopic):
    descriptor = descriptor.strip()
    descriptor = _AMBIGUITY_SUFFIX.sub("", descriptor).strip()
    if not descriptor or descriptor.startswith("#"):
        return None
    lowered = descriptor.lower()
    if lowered.startswith("info:") or lowered.startswith("obs:"):
        return None
    if lowered.startswith(("glycan:", "gno:", "gnome:", "xlmod:", "xmod:")):
        raise ValueError("GNO/GNOme and XL-MOD modifications are not yet supported")
    if _NUMBER.fullmatch(descriptor):
        return float(descriptor)
    if lowered.startswith("formula:"):
        return _formula_mass(descriptor.split(":", 1)[1], monoisotopic)

    prefix, separator, value = descriptor.partition(":")
    normalized_prefix = prefix.upper()
    if separator and normalized_prefix in {"U", "UNIMOD"}:
        if value.isdigit():
            return _lookup("UNIMOD:" + value, site, monoisotopic)
        return _lookup(value, site, monoisotopic, "UNIMOD")
    if separator and normalized_prefix in {"M", "MOD", "PSI-MOD"}:
        if value.isdigit():
            return _lookup("MOD:" + value.zfill(5), site, monoisotopic)
        return _lookup(value, site, monoisotopic, "MOD")
    if separator and normalized_prefix in {"R", "RESID"}:
        if value.upper().startswith("AA"):
            return _lookup("RESID:" + value.upper(), site, monoisotopic)
        return _lookup(value, site, monoisotopic, "RESID")
    return _lookup(descriptor, site, monoisotopic)


def _tag_mass(tag, site, monoisotopic):
    if tag.strip().startswith("#"):
        return 0.0
    values = []
    errors = []
    for descriptor in tag.split("|"):
        try:
            value = _descriptor_mass(descriptor, site, monoisotopic)
        except ValueError as exception:
            errors.append(exception)
            continue
        if value is not None:
            values.append(value)
    if not values:
        if errors:
            raise errors[0]
        raise ValueError("Modification tag has no mass-bearing descriptor: {!r}".format(tag))
    if not all(math.isclose(values[0], value, abs_tol=0.01) for value in values[1:]):
        raise ValueError("Conflicting masses in modification tag: {!r}".format(tag))
    return values[0]


def calc_proforma_mass(sequence, monoisotopic=True, ion_type="H2O"):
    """Calculate a neutral mass for a supported ProForma protein string."""
    if __package__:
        from . import mass as base_mass
    else:
        import mass as base_mass

    plain, modifications, x_gaps, terminal_modification = _parse(sequence)
    if terminal_modification and str(ion_type).upper() != "H2O":
        raise ValueError("Explicit terminal modifications cannot be combined with fragment ion_type")

    masses = (
        base_mass.aa_masses_monoisotopic.copy()
        if monoisotopic else base_mass.aa_masses.copy()
    )
    masses.update({
        "J": masses["I"],
        "B": (masses["D"] + masses["N"]) / 2,
        "Z": (masses["E"] + masses["Q"]) / 2,
        "O": 237.147727 if monoisotopic else 237.29816,
        "U": 150.953634 if monoisotopic else 150.0379,
        "X": 0.0,
    })
    for tags in x_gaps:
        if not tags or not any(
            _NUMBER.fullmatch(_AMBIGUITY_SUFFIX.sub("", part).strip())
            for tag in tags for part in tag.split("|")
        ):
            raise ValueError("X requires an explicit signed mass gap, for example X[+367.0537]")

    total = sum(masses[residue] for residue in plain)
    total += base_mass.get_pep_ion_mass_shift(ion_type, monoisotopic=monoisotopic)
    for tag, site in modifications:
        total += _tag_mass(tag, site, monoisotopic)
    return float(total)
