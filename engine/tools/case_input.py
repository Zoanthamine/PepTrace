"""Shared input checks and canonical linear peptide representation."""
from knowledge_source import components, attachment_table
from mass_catalog import group_formula
from mass_ledger import MONO, RESIDUE

C_TERM_DELTAS = components()['c_terminal_deltas']

SULFUR_BOND_DELTAS = components()['sulfur_bond_deltas']


def check_keys(value, allowed, label):
    if not isinstance(value, dict):
        raise ValueError(f"{label} requires an object")
    unknown = value.keys() - allowed
    if unknown:
        raise ValueError(f"Unknown {label} fields: {sorted(unknown)}")


def string_set(value, label):
    if not isinstance(value, list) or any(not isinstance(item, str) or not item.strip() for item in value):
        raise ValueError(f"{label} requires a list of nonempty strings")
    return set(value)


def read_sequence(sequence):
    residues = list(sequence) if isinstance(sequence, str) else sequence
    if not isinstance(residues, list) or not residues or any(
        not isinstance(r, str) or r not in RESIDUE or r == "O" for r in residues
    ):
        raise ValueError("Use uppercase one-letter residues or a list including named residues such as Orn")
    return list(residues)


# Curated attachment scopes for protection groups on the supported linear model.
# Aggregate or custom elemental changes remain explicit formula bookkeeping.
ATTACHMENT = attachment_table()


def normalize_peptide(record):
    check_keys(record, {"sequence", "nterm", "cterm", "modifications", "topology", "fragment_model", "bonds", "alpha_stereochemistry", "monomer_state"}, "peptide")
    residues = read_sequence(record["sequence"])
    nterm, cterm = record["nterm"], record["cterm"]
    if nterm not in ("free", "acetyl") or not isinstance(cterm, str) or cterm not in C_TERM_DELTAS:
        raise ValueError(f"Supported termini: free/acetyl and {', '.join(C_TERM_DELTAS)}")
    if record.get("topology", "linear") not in {"linear", "crosslinked"}:
        raise ValueError("Supported input is one backbone chain, optionally with explicit sulfur bridges")
    if record.get("fragment_model", "linear_local_modifications") != "linear_local_modifications":
        raise ValueError("Unsupported fragment model")
    mods = record.get("modifications", [])
    if not isinstance(mods, list):
        raise ValueError("modifications requires a list")
    normalized = []
    occupied_mods = set()
    for mod in mods:
        check_keys(mod, {"site", "group", "delta_formula", "name"}, "modification")
        if "name" in mod and (not isinstance(mod["name"], str) or not mod["name"].strip()):
            raise ValueError("Modification name must be a nonempty display label")
        site = mod["site"]
        if (type(site) is not int or not 1 <= site <= len(residues)) and site not in ("N", "C"):
            raise ValueError(f"Invalid modification site: {site}")
        if ("group" in mod) == ("delta_formula" in mod):
            raise ValueError("Specify exactly one of group or delta_formula for each modification")
        if (site == "N" and nterm != "free") or (site == "C" and cterm != "acid"):
            raise ValueError("Terminal modification overlaps an explicit terminal state")
        group = mod.get("group")
        if "group" in mod:
            if not isinstance(group, str):
                raise ValueError("group must be a catalog string")
            group_formula(group)
            if group == "NalphaAzide" and residues[0] == "P":
                raise ValueError("NalphaAzide requires a primary N-terminal amine; Pro is secondary")
            residue = f"{site}-term" if isinstance(site, str) else residues[site - 1]
            allowed = ATTACHMENT.get(group)
            if allowed is not None and residue not in allowed:
                raise ValueError(f"{group} is incompatible with site {site} in the linear model")
        else:
            formula = mod["delta_formula"]
            if not isinstance(formula, dict) or not formula or any(
                element not in MONO or type(count) is not int for element, count in formula.items()
            ):
                raise ValueError("delta_formula requires an element-to-integer mapping")
        if site in occupied_mods:
            raise ValueError(f"Overlapping modifications at site {site}; specify one net change or an explicit structure")
        occupied_mods.add(site)
        if group == "C-amide":
            if site != "C":
                raise ValueError("C-amide requires a free acid C terminus")
            cterm = "amide"
        elif group == "Ac" and site == "N":
            nterm = "acetyl"
        else:
            normalized.append(dict(mod))
    bonds = record.get("bonds", [])
    if not isinstance(bonds, list):
        raise ValueError("bonds requires a list")
    bonded = set()
    canonical_bonds = []
    for bond in bonds:
        check_keys(bond, {"kind", "sites", "oxidized_site"}, "bond")
        kind, sites = bond["kind"], bond["sites"]
        if kind not in SULFUR_BOND_DELTAS:
            raise ValueError(f"Supported sulfur bonds: {', '.join(SULFUR_BOND_DELTAS)}")
        if not isinstance(sites, list) or len(sites) != 2 or any(
            type(i) is not int or not 1 <= i <= len(residues) or residues[i-1] != "C" for i in sites
        ) or sites[0] == sites[1]:
            raise ValueError("A sulfur bridge requires two distinct Cys sites")
        if set(sites) & (bonded | occupied_mods):
            raise ValueError("Sulfur bond overlaps another bridge or side-chain modification")
        bonded.update(sites)
        if kind in {"thiosulfinate", "thiosulfonate"}:
            if bond.get("oxidized_site") not in sites:
                raise ValueError("An oxidized sulfur bridge requires its explicit oxidized_site")
        elif "oxidized_site" in bond:
            raise ValueError("This sulfur bridge cannot carry oxidized_site")
        canonical_bonds.append(dict(bond, sites=sorted(sites)))
    result = {"sequence": residues, "nterm": nterm, "cterm": cterm, "modifications": normalized}
    complex_groups={m.get("group") for m in normalized}&{"NCA","NChloroformylAcidChloride"}
    if complex_groups:
        if len(residues)!=1 or cterm!="acid" or any(m["site"]=="C" for m in normalized) or residues[0] in {"beta-Ala","Aoa"}:
            raise ValueError("NCA/chloroformyl-chloride state requires one alpha amino acid and no separate C-terminal state")
        state={"NCA":"amino acid NCA","NChloroformylAcidChloride":"N-chloroformyl amino acid chloride"}[next(iter(complex_groups))]
        if "monomer_state" in record and record["monomer_state"]!=state:
            raise ValueError("monomer_state disagrees with the explicit combined terminal group")
        result["monomer_state"]=state
    elif "monomer_state" in record:
        raise ValueError("monomer_state requires an explicit NCA/chloroformyl-chloride group")
    stereo = record.get("alpha_stereochemistry", {})
    if not isinstance(stereo, dict) or any(
        not isinstance(key, str) or not key.isdigit() or str(int(key)) != key
        or not 1 <= int(key) <= len(residues) or value not in {"L", "D"}
        or residues[int(key)-1] in {"G", "Aib", "beta-Ala", "Aoa", "Dha", "MeDha"}
        for key, value in stereo.items()
    ):
        raise ValueError("alpha_stereochemistry requires 1-based string sites and L/D at chiral alpha centers")
    if stereo:
        result["alpha_stereochemistry"] = dict(stereo)
    if canonical_bonds:
        result.update(bonds=sorted(canonical_bonds, key=lambda b: tuple(b["sites"])), topology="crosslinked")
    elif record.get("topology") == "crosslinked":
        raise ValueError("crosslinked input requires explicit bonds")
    return result


def evaluate_requirements(all_of, any_of, present, absent):
    """Return support, missing requirements and conflicts for one fact family."""
    support = set(all_of) & present
    missing = set(all_of) - present - absent
    conflicts = set(all_of) & absent
    if any_of:
        alternatives = set(any_of)
        label = " or ".join(sorted(alternatives))
        if alternatives & present:
            support.update(alternatives & present)
        elif alternatives <= absent:
            conflicts.add(label)
        else:
            missing.add(label)
    return sorted(support), sorted(missing), sorted(conflicts)


def add_starting_group_facts(starting, target, context):
    """Derive retention prerequisites only when the residue identities coincide."""
    if starting is None or target is None or starting["sequence"]!=target["sequence"]:
        return
    present, absent=context
    checks=[("n_fmoc_precursor:N","N","Fmoc")]
    checks += [(f"trp_formyl_protected:{i}",i,"Formyl") for i,r in enumerate(starting["sequence"],1) if r=="W"]
    checks += [(f"met_oxidized_precursor:{i}",i,"Oxidation") for i,r in enumerate(starting["sequence"],1) if r=="M"]
    for fact,site,group in checks:
        exists=any(m["site"]==site and m.get("group")==group for m in starting["modifications"])
        yes,no=(present,absent) if exists else (absent,present)
        if fact in no["features"]:
            raise ValueError(f"Starting material contradicts context feature {fact}")
        yes["features"].add(fact)
