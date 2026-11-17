"""Find source-backed leads for a reported neutral mass difference."""
import math
import re

from mass_ledger import RESIDUE, formula_mass, has_neutral_lookup_delta
from case_input import check_keys, evaluate_requirements, read_sequence, string_set
from knowledge_source import execution_adapters

RELATIVE_KINDS = {"nominal_delta", "delta_mono", "delta_approx"}
BOOK_NOMINAL = re.compile(r"^\s*([+\-−–]?\d+)\s*amu\b", re.I)
TARGET_REFERENCE = "intended_product"

# Reviewed applicability for relative search. Untyped entries remain unfiltered leads.
# The keys are existing library rules. Values are required sets of context facts.
PROFILES = {rid:{k:set(v) for k,v in profile.items()} for rid,profile in execution_adapters()['relative_profiles'].items()}
STAGE_ALTERNATIVES = {rid:set(v) for rid,v in execution_adapters()['stage_alternatives'].items()}


def read_context(case, tags, target_sequence=None):
    context = case.get("context", {})
    if not isinstance(context, dict):
        raise ValueError("context requires an object")
    unknown = context.keys() - {"sequence", "residues_present", "features_present", "features_absent", "stage_tags_absent"}
    if unknown:
        raise ValueError(f"Unknown context keys: {sorted(unknown)}")
    sequence = read_sequence(context["sequence"]) if "sequence" in context else target_sequence
    if target_sequence is not None and sequence != target_sequence:
        raise ValueError("context.sequence contradicts intended_product.sequence")
    named = string_set(context.get("residues_present", []), "residues_present")
    residues = set(sequence) if sequence is not None else named
    if not residues <= RESIDUE.keys() - {"O"} or not named <= residues:
        raise ValueError("Unknown or inconsistent residue facts")
    present = {"stage_tags": set(tags), "residues": residues,
               "features": string_set(context.get("features_present", []), "features_present")}
    absent_tags = string_set(case["reaction"].get("stage_tags_absent", []), "reaction.stage_tags_absent")
    absent = {"stage_tags": absent_tags | string_set(context.get("stage_tags_absent", []), "context.stage_tags_absent"),
              "residues": RESIDUE.keys() - residues if sequence is not None else set(),
              "features": string_set(context.get("features_absent", []), "features_absent")}
    if any(present[key] & absent[key] for key in present):
        raise ValueError("Context contains contradictory facts")
    return present, absent


def stage_fit(rule_id, context):
    present, absent = context
    return evaluate_requirements(PROFILES.get(rule_id, {}).get("stage_tags", ()),
                                 STAGE_ALTERNATIVES.get(rule_id, ()),
                                 present["stage_tags"], absent["stage_tags"])


def applicability(rule_id, context, profile=None):
    if profile is None and rule_id not in PROFILES:
        return "unreviewed_context", [], [], []
    profile=profile if profile is not None else PROFILES[rule_id]
    present, absent = context
    supported, missing, conflicting = evaluate_requirements(profile.get("stage_tags",[]),
        profile.get("stage_any",STAGE_ALTERNATIVES.get(rule_id,[])),present["stage_tags"],absent["stage_tags"])
    for key, facts in profile.items():
        if key in {"residues","features"}:
            support, unknown, conflict = evaluate_requirements(facts, (), present[key], absent[key])
            supported.extend(support)
            missing.extend(unknown)
            conflicting.extend(conflict)
    support,unknown,conflict=evaluate_requirements([],profile.get("any_residues",[]),present["residues"],absent["residues"])
    supported.extend(support)
    missing.extend(unknown)
    conflicting.extend(conflict)
    if conflicting:
        return "context_conflict", supported, missing, conflicting
    return ("context_partial" if missing else "context_supported"), supported, missing, conflicting


def validate_relative(item):
    value = item["value"]
    reference = item["reference"]
    nominal = item["kind"] == "nominal_delta"
    if reference not in (TARGET_REFERENCE, "starting_material") or (nominal and reference != TARGET_REFERENCE):
        raise ValueError("Nominal source lookup requires a difference versus the intended product")
    if nominal:
        if type(value) is not int or {"tolerance", "tolerance_unit", "charge", "ions"} & item.keys():
            raise ValueError("nominal_delta requires an integer and no exact-mass or ion fields")
    elif (item.get("tolerance_unit") != "Da" or not all(
            type(n) in (int, float) and math.isfinite(n) for n in (value, item.get("tolerance"))
            ) or item["tolerance"] <= 0 or {"charge", "ions"} & item.keys()):
        raise ValueError(f'{item["kind"]} requires a finite neutral difference and positive Da tolerance')
    check_keys(item, {"kind", "value", "reference", "tolerance", "tolerance_unit"}, "relative observation")


def search_relative(case, rules, context, appendix=None):
    for item in case["observations"]:
        validate_relative(item)
        if item["reference"] != TARGET_REFERENCE:
            yield {"input": item, "matches": [], "reference_only_delta_leads": [],
                   "appendix_only_leads": [], "status": "no_relative_lead"}
            continue
        value = item["value"]
        nominal = item["kind"] == "nominal_delta"
        approximate = item["kind"] == "delta_approx"

        leads = []
        for rid, rule in rules.items():
            if rule.get("automatic_mass_event") is False and not rule.get("nominal_reference_lookup_allowed"):
                continue
            status = rule["mass_status"]
            if status == "not_product_mass_rule":
                continue
            mass = (0.0 if status == "isomeric" else formula_mass(rule["delta_formula"])) if has_neutral_lookup_delta(rule) else None
            label = BOOK_NOMINAL.match(rule.get("nominal_mass_shift_book") or "")
            if nominal:
                if label:
                    matched = int(label[1].replace("−", "-").replace("–", "-")) == value
                    basis = "book_reported_integer"
                else:
                    matched = mass is not None and abs(mass - value) < .5
                    basis = "rounded_theoretical_neutral_delta"
            elif mass is not None and abs(value - mass) <= item["tolerance"]:
                matched = True
                basis = "theoretical_neutral_delta"
            elif approximate and label is not None:
                matched = abs(value - int(label[1].replace("−", "-").replace("–", "-"))) <= item["tolerance"]
                basis = "book_reported_integer_in_window"
            else:
                matched = mass is not None and abs(value - mass) <= item["tolerance"]
                basis = "theoretical_neutral_delta"
            if not matched:
                continue
            fit, support, missing, conflicts = applicability(rid, context, rule.get("relative_profile"))
            leads.append({"rule_id": rid, "name_cn": rule["name_cn"], "origin": rule["origin"],
                          "source": rule.get("source", rule.get("source_basis")),
                          "mass_basis": basis, "mass_status": status, "theoretical_delta_mono": mass,
                          "delta_formula": rule.get("delta_formula") if mass is not None else None,
                          "book_nominal_label": rule.get("nominal_mass_shift_book"),
                          "error_da": value - mass if basis == "theoretical_neutral_delta" else None,
                          "reference": rule.get("mass_reference"), "context_status": fit,
                          "support": support, "missing": missing, "conflicts": conflicts,
                          "verification": rule.get("diagnostic_evidence", []),
                          "structure_status": "unconfirmed_hypothesis" if mass is None else "rule_based_candidate"})
        priority = ["context_supported", "context_partial", "unreviewed_context", "context_conflict"]
        leads.sort(key=lambda lead: (priority.index(lead["context_status"]), -len(lead["support"]),
                                    lead["mass_basis"] == "rounded_theoretical_neutral_delta", lead["rule_id"]))
        appendix_leads = []
        for source, ledger in appendix or ():
            mass = ledger["exact_neutral_delta"]
            number = source["nominal_average_delta_book"]
            matched = (number == value if nominal else
                       mass is not None and abs(value - mass) <= item["tolerance"] or
                       approximate and number is not None and abs(value - number) <= item["tolerance"])
            if matched:
                appendix_leads.append({
                    "appendix_id": source["appendix_row_id"], "name": source["modification_text"],
                    "print_page": source["print_page"], "pdf_page": source["pdf_page"],
                    "book_nominal": source["nominal_average_delta_book"],
                    "theoretical_delta_mono": mass, "reference_note": ledger["note"],
                    **{key: ledger[key] for key in ("mass_reference_kind", "intrinsic_charge",
                                                   "elemental_mass_delta", "same_protonated_ion_delta")
                       if key in ledger},
                    "evidence": "appendix_reference_only",
                })
        yield {"input": item, "matches": [], "reference_only_delta_leads": leads,
               "appendix_only_leads": appendix_leads,
               "status": "relative_leads_only" if leads or appendix_leads else "no_relative_lead"}
