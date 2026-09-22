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


def _parse_core(text, offset=0):
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
            position = offset + len(sequence) - 1
            index += 1
            tags = []
            while index < len(text) and text[index] == "[":
                end = _matching(text, index, "[", "]")
                tag = text[index + 1:end]
                modifications.append((tag, residue, frozenset({position})))
                tags.append(tag)
                index = end + 1
            if residue == "X":
                x_gaps.append(tags)
            continue
        if character == "(":
            end = _matching(text, index, "(", ")")
            region_start = offset + len(sequence)
            region_sequence, region_mods, region_x, region_terminal = _parse_core(
                text[index + 1:end], region_start
            )
            sequence.extend(region_sequence)
            modifications.extend(region_mods)
            x_gaps.extend(region_x)
            has_terminal_modification |= region_terminal
            index = end + 1
            while index < len(text) and text[index] == "[":
                tag_end = _matching(text, index, "[", "]")
                positions = frozenset(
                    range(region_start, region_start + len(region_sequence))
                )
                modifications.append((text[index + 1:tag_end], None, positions))
                index = tag_end + 1
            continue
        if character == "-" and index + 1 < len(text) and text[index + 1] == "[":
            has_terminal_modification = True
            index += 1
            while index < len(text) and text[index] == "[":
                end = _matching(text, index, "[", "]")
                modifications.append((text[index + 1:end], "C-term", "C-term"))
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
        modifications.append((text[1:end], None, None))
        text = text[end + 1:]

    while text.startswith("["):
        end = _matching(text, 0, "[", "]")
        remainder = text[end + 1:]
        occurrence = re.match(r"\^(\d+)\?", remainder)
        if remainder.startswith("?") or occurrence:
            count = int(occurrence.group(1)) if occurrence else 1
            if count < 1:
                raise ValueError("Unlocalized modification count must be positive")
            modifications.extend((text[1:end], None, None) for _ in range(count))
            text = remainder[occurrence.end():] if occurrence else remainder[1:]
        elif remainder.startswith("-"):
            modifications.append((text[1:end], "N-term", "N-term"))
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
                positions = [
                    index for index, residue in enumerate(base_sequence)
                    if residue == selector.upper()
                ]
                modifications.extend(
                    (tag, selector.upper(), frozenset({position}))
                    for position in positions
                )
            elif normalized in {"n-term", "c-term"}:
                site = normalized[0].upper() + normalized[1:]
                modifications.append((tag, site, site))
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


def _fragment_mass_options(sequence, modifications, ion_type, length, monoisotopic):
    """Return the distinct masses possible for one backbone fragment."""
    if __package__:
        from . import mass as base_mass
    else:
        import mass as base_mass

    masses = (
        base_mass.aa_masses_monoisotopic
        if monoisotopic else base_mass.aa_masses
    )
    residue_options = {
        **{residue: (mass,) for residue, mass in masses.items()},
        "J": (masses["I"],),
        "B": (masses["D"], masses["N"]),
        "Z": (masses["E"], masses["Q"]),
        "O": (237.147727 if monoisotopic else 237.29816,),
        "U": (150.953634 if monoisotopic else 150.0379,),
        "X": (0.0,),
    }
    size = len(sequence)
    positions = (
        frozenset(range(length))
        if ion_type[0] in {"a", "b", "c"}
        else frozenset(range(size - length, size))
    )
    options = [base_mass.get_pep_ion_mass_shift(
        ion_type, monoisotopic=monoisotopic
    )]
    for position in sorted(positions):
        options = sorted({
            round(total + residue_mass, 12)
            for total in options
            for residue_mass in residue_options[sequence[position]]
        })

    groups = {}
    for tag, _, candidate_positions in modifications:
        match = _AMBIGUITY_SUFFIX.search(tag)
        if match and isinstance(candidate_positions, frozenset):
            groups.setdefault(match.group(1), set()).update(candidate_positions)

    modification_options = []
    for tag, site, candidate_positions in modifications:
        if tag.strip().startswith("#"):
            continue
        is_unlocalized = candidate_positions is None
        match = _AMBIGUITY_SUFFIX.search(tag)
        if match:
            candidate_positions = frozenset(groups[match.group(1)])
        elif candidate_positions is None:
            candidate_positions = frozenset(range(size))
        modification_options.append(
            (_tag_mass(tag, site, monoisotopic), candidate_positions, is_unlocalized)
        )

    # Repeated unlocalized modifications share a site pool and cannot occupy
    # more candidate sites than exist on either side of a cleavage.
    counted = {}
    for mass, candidates, is_unlocalized in modification_options:
        key = (mass, candidates, is_unlocalized)
        counted[key] = counted.get(key, 0) + 1
    for (modification_mass, candidates, is_unlocalized), count in counted.items():
        if candidates == "N-term":
            counts = (count,) if ion_type[0] in {"a", "b", "c"} else (0,)
        elif candidates == "C-term":
            counts = (count,) if ion_type[0] in {"x", "y", "z"} else (0,)
        else:
            inside = len(positions & candidates)
            outside = len(candidates) - inside
            if is_unlocalized:
                minimum = max(0, count - outside)
                maximum = min(count, inside)
                counts = range(minimum, maximum + 1)
            elif not inside:
                counts = (0,)
            elif not outside:
                counts = (count,)
            else:
                counts = range(count + 1)
        options = sorted({
            round(total + occurrence_count * modification_mass, 12)
            for total in options
            for occurrence_count in counts
        })

    return tuple(sorted({round(value, 12) for value in options}))


def _calc_proforma_fragments(
    sequence,
    ion_types=("b", "y"),
    monoisotopic=True,
    ambiguous_rule="reject",
):
    """Calculate backbone fragments for a supported ProForma sequence."""
    if __package__:
        from . import mass as base_mass
    else:
        import mass as base_mass

    plain, modifications, _, _ = _parse(sequence)
    # Validate mass gaps and every modification even if all resulting
    # fragments are rejected as ambiguous.
    calc_proforma_mass(sequence, monoisotopic=monoisotopic)

    fragments = {}
    for ion_type in ion_types:
        for length in range(1, len(plain)):
            options = _fragment_mass_options(
                plain, modifications, ion_type, length, monoisotopic
            )
            fragment_name = base_mass._pep_fragment_name(ion_type, length)
            if len(options) == 1:
                fragments[fragment_name] = options[0]
            elif ambiguous_rule == "both":
                for index, mass in enumerate(options, 1):
                    fragments[f"{fragment_name}#{index}"] = mass
    return fragments


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
    for tag, site, _ in modifications:
        total += _tag_mass(tag, site, monoisotopic)
    return float(total)


if __name__ == "__main__":

    ca_seq = "S[Acetylation]HHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEATSRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLVHWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNVLDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRGFPK"

    experimental_mass_text = """374.208
    420.187
    587.331
    588.339
    605.259
    606.266
    663.288
    702.418
    740.987
    816.463
    883.373
    970.559
    971.556
    1011.468
    1084.602
    1104.513
    1124.506
    1147.944
    1148.527
    1212.697
    1218.556
    1262.57
    1325.781
    1416.644
    1454.769
    1545.687
    1550.893
    1638.733
    1682.746
    1766.747
    1868.826
    1875.083
    1961.87
    2005.884
    2061.163
    2133.979
    2175.206
    2241.105
    2246.243
    2249.006
    2257.079
    2257.581
    2352.058
    2493.127
    2549.311
    2606.212
    2633.235
    2677.248
    2716.536
    2791.292
    2844.706
    2845.579
    2847.434
    2847.936
    2848.312
    2850.319
    2853.938
    2886.481
    2977.356
    2981.132
    3070.025
    3070.527
    3071.674
    3087.693
    3088.963
    3089.443
    3128.696
    3135.454
    3144.715
    3217.501
    3226.125
    3257.74
    3261.516
    3312.061
    3312.395
    3328.777
    3344.795
    3398.202
    3398.704
    3442.82
    3445.6
    3457.834
    3460.845
    3462.48
    3483.734
    3484.236
    3544.67
    3589.89
    3614.679
    3624.105
    3659.696
    3701.161
    3703.932
    3714.014
    3714.516
    3726.336
    3726.837
    3728.767
    3753.981
    3772.781
    3887.807
    3892.287
    3899.868
    3899.917
    3907.927
    3918.063
    3921.086
    3988.856
    4074.169
    4087.167
    4090.163
    4116.952
    4173.695
    4190.229
    4204.877
    4221.238
    4237.266
    4386.125
    4463.421
    4472.202
    4514.14
    4514.186
    4565.471
    4585.196
    4721.515
    4725.26
    4741.831
    4753.287
    4797.274
    4797.302
    4798.33
    4810.313
    4888.579
    4910.358
    4910.386
    4936.608
    4939.648
    4940.658
    5023.639
    5062.596
    5120.158
    5123.726
    5135.534
    5138.741
    5147.695
    5153.686
    5200.756
    5209.736
    5241.763
    5248.618
    5319.656
    5337.834
    5338.639
    5428.855
    5432.74
    5486.789
    5490.843
    5531.808
    5548.918
    5565.96
    5589.949
    5672.025
    5677.013
    5678.045
    5694.872
    5706.875
    5709.871
    5729.292
    5751.892
    5823.321
    5824.913
    5875.131
    5880.937
    5905.185
    5950.972
    6002.252
    6053.02
    6107.24
    6107.745
    6139.047
    6142.984
    6164.755
    6194.282
    6235.516
    6296.153
    6299.055
    6452.253
    6475.484
    6514.486
    6604.522
    6614.031
    6682.364
    6708.14
    6717.61
    6752.396
    6790.608
    6795.402
    6831.7
    6909.445
    6967.472
    7089.259
    7104.532
    7125.845
    7190.562
    7225.884
    7227.916
    7232.034
    7257.846
    7294.619
    7297.578
    7317.396
    7321.963
    7331.544
    7337.627
    7339.976
    7340.642
    7340.84
    7344.519
    7347.54
    7352.541
    7408.66
    7427.006
    7452.674
    7485.05
    7494.938
    7508.698
    7550.736
    7551.748
    7610.03
    7627.218
    7632.028
    7663.761
    7680.786
    7745.18
    7769.767
    7799.837
    7800.122
    7843.851
    7846.218
    7941.875
    7958.877
    8029.89
    8032.29
    8073.904
    8078.94
    8160.937
    8171.516
    8194.348
    8210.385
    8258.36
    8288.994
    8293.374
    8309.367
    8318.641
    8405.034
    8422.455
    8521.525
    8532.118
    8574.517
    8603.154
    8635.53
    8635.573
    8639.045
    8640.531
    8702.22
    8800.686
    8815.307
    8817.722
    8829.7
    8845.699
    8861.746
    8943.408
    8958.784
    9045.851
    9058.431
    9271.918
    9297.916
    9313.847
    9314.913
    9426.629
    9462
    9463.06
    9483.655
    9519.042
    9577.052
    9584.699
    9690.124
    9747.768
    9820.149
    9837.195
    9892.175
    9893.201
    9903.867
    9913.987
    9918.188
    9936.026
    9936.202
    9952.224
    9954.154
    10021.24
    10037.252
    10079.271
    10107.285
    10116.018
    10123.269
    10139.273
    10238.827
    10239.327
    10242.071
    10244.086
    10251.371
    10308.378
    10391.148
    10437.498
    10494.538
    10512.188
    10520.523
    10528.207
    10537.521
    10554.562
    10666.644
    10768.322
    10779.722
    10795.736
    10812.332
    10822.748
    10849.745
    10865.72
    10881.756
    10980.756
    10980.765
    10999.416
    10999.461
    11067.836
    11069.856
    11165.906
    11180.927
    11181.905
    11229.497
    11280.941
    11296.958
    11344.529
    11459.558
    11461.567
    11464.94
    11544.606
    11554.903
    11578.527
    11588.616
    11620.18
    11644.62
    11695.866
    11731.665
    11731.668
    11860.717
    11862.348
    11862.348
    11873.382
    11932.346
    11997.745
    11999.995
    12100.447
    12127.441
    12128.456
    12143.431
    12182.873
    12197.887
    12215.507
    12312.911
    12330.529
    12345.541
    12348.565
    12373.858
    12386.525
    12424.995
    12425.001
    12466.043
    12469.017
    12485.6
    12576.742
    12596.1
    12596.104
    12614.701
    12683.791
    12727.795
    12829.259
    12839.828
    12844.257
    12874.862
    12877.689
    12882.445
    12888.265
    12959.307
    12972.897
    12987.938
    12989.95
    13026.976
    13029.339
    13030.952
    13131.044
    13159.386
    13228.044
    13231.733
    13256.113
    13301.146
    13410.591
    13411.592
    13413.214
    13470.225
    13533.396
    13543.164
    13566.25
    13584.228
    13661.388
    13701.633
    13750.329
    13760.803
    13766.373
    13792.354
    13810.351
    13938.428
    13945.82
    13954.562
    13992.484
    14010.479
    14058.869
    14079.491
    14159.904
    14180.559
    14238.593
    14271.98
    14288.008
    14300.541
    14353.585
    14385.643
    14451.07
    14451.077
    14455.782
    14457.689
    14466.077
    14482.648
    14483.666
    14484.644
    14500.681
    14504.367
    14504.874
    14509.096
    14514.674
    14516.674
    14541.691
    14557.707
    14587.947
    14607.902
    14625.127
    14660.747
    14676.758
    14702.753
    14720.766
    14735.762
    14818.994
    14827.212
    14838.347
    14848.85
    14862.835
    14943.16
    15000.296
    15003.304
    15019.93
    15064.966
    15185.404
    15189.467
    15198.398
    15202.46
    15206.024
    15250.042
    15263.407
    15265.027
    15265.378
    15266.059
    15270.542
    15294.422
    15310.43
    15314.491
    15424.514
    15529.092
    15538.528
    15599.204
    15697.495
    15707.625
    15711.643
    15736.297
    15850.393
    15861.672
    15879.742
    15881.697
    15908.138
    15957.699
    15977.399
    15978.811
    15998.423
    16032.814
    16033.817
    16036.83
    16048.435
    16050.475
    16061.842
    16120.496
    16135.495
    16186.36
    16284.575
    16410.645
    16494.723
    16539.758
    16653.846
    16678.867
    16678.912
    16695.861
    16710.711
    16766.867
    16766.871
    16795.836
    16810.844
    16863.299
    16864.311
    17012.006
    17063.229
    17148.047
    17364.139
    17403.447
    17549.188
    17665.229
    18144.008
    18196.418
    18480.533
    18620.605
    18687.611
    18766.707
    18893.729
    18965.295
    19072.543
    19160.045
    19198.357
    19619.334
    19693.82
    19949.369
    19950.32
    20307.54
    20476.635
    20935.838
    21049.867
    21072.861
    21336.203
    21651.549
    21885.939
    22212.355
    22227.01
    22379.918
    22390.701
    22418.822
    22688.047
    22790.852
    22872.795
    23393.357
    23722.426
    24429.949
    24706.971
    24729.447
    25449.395
    25451.502
    26509.49
    26534.086
    26597.18
    27802.301
    27949.738
    28131.293
    28266.793
    28301.156
    28308.416
    28765.414
    28945.748
    28949.719
    28952.771
    28965.74
    28968.617
    28969.766
    28991.754
    29008.748
    29009.773
    29195.82
    29216.879
    29222.727
    29231.93
    29245.232
    29284.625
    29527.086
    29535.311
    29715.912
    30398.807
    30837.537
    31416.27
    31695.678
    31955.619
    32074.863
    32957.918
    43498.945
    """

    experimental_masses = [
        float(value) for value in experimental_mass_text.split()
    ]
    predicted_fragments = _calc_proforma_fragments(
        ca_seq, ion_types=("c", "z'"), monoisotopic=True
    )
    tolerance_ppm = 20.0
    match_candidates = sorted(
        (
            abs(experimental - predicted) / predicted * 1e6,
            ion,
            predicted,
            experimental,
        )
        for ion, predicted in predicted_fragments.items()
        for experimental in experimental_masses
        if abs(experimental - predicted) / predicted * 1e6 <= tolerance_ppm
    )

    # Select the lowest-error one-to-one assignments so neither an ion nor an
    # experimental feature can inflate the sequence coverage.
    matches = []
    matched_ions = set()
    matched_features = set()
    for _, ion, predicted, experimental in match_candidates:
        if ion in matched_ions or experimental in matched_features:
            continue
        matched_ions.add(ion)
        matched_features.add(experimental)
        error_ppm = (experimental - predicted) / predicted * 1e6
        matches.append((ion, predicted, experimental, error_ppm))

    sequence_length = len(strip_proforma(ca_seq))
    cleavage_sites = {
        int(ion.lstrip("abcxyz'").partition("#")[0])
        if ion.startswith("c")
        else sequence_length - int(ion.lstrip("abcxyz'").partition("#")[0])
        for ion, _, _, _ in matches
    }
    coverage = len(cleavage_sites) / (sequence_length - 1) * 100

    print("Ion\tPredicted mass\tExperimental mass\tError (ppm)")
    for ion, predicted, experimental, error_ppm in sorted(
        matches,
        key=lambda match: (
            match[0][0],
            int(match[0].lstrip("abcxyz'").partition("#")[0]),
        ),
    ):
        print(f"{ion}\t{predicted:.6f}\t{experimental:.3f}\t{error_ppm:+.3f}")
    print(
        f"Matched {len(matches)} fragments at {tolerance_ppm:g} ppm; "
        f"{len(cleavage_sites)}/{sequence_length - 1} cleavage sites; "
        f"sequence coverage {coverage:.2f}%"
    )
