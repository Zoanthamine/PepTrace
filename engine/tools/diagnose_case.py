"""Evidence-aware peptide diagnosis. CLI input and output are UTF-8 JSON."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

from mass_ledger import MONO, PROTON, RESIDUE, formula_mass, has_neutral_lookup_delta, parse_formula
from mass_catalog import group_formula
from case_input import C_TERM_DELTAS, SULFUR_BOND_DELTAS, add_starting_group_facts, check_keys, normalize_peptide, string_set
from linear_rules import instantiate, read_rule_pack, sequence_rule_exclusions
from sulfur_bond_rules import instantiate_bond_rules
from carbonyl_components import normalize_carbonyl_donors
from relative_case import RELATIVE_KINDS, applicability, read_context, search_relative, stage_fit
from knowledge_source import execution_adapters, load_runtime
from observation_model import parse_coordinate, formula_mass as observation_formula_mass

BOOK = Path(__file__).resolve().parents[1]
ELECTRON = 0.000548579909
IONS = {"H": (1, PROTON), "-H": (-1, -PROTON),
        "Na": (1, MONO["Na"] - ELECTRON), "K": (1, MONO["K"] - ELECTRON)}
# Explicit activation criteria. Unreviewed book rules remain reference leads.
SITE_RULES = execution_adapters()['site_rules']
RETENTION_RULES = execution_adapters()['retention_rules']
REVIEWED_IDS = set(execution_adapters()['reviewed_ids'])


def load_rules():
    bundle, version = load_runtime()
    by_id={r['rule_id']:dict(r['legacy_record'],origin=r['legacy_record'].get('origin',r['origin']))
           for r in bundle['rules'] if r['legacy_record'] is not None}
    for record in bundle['rules']:
        if record['legacy_record'] is not None and record.get('effective_source'):
            by_id[record['rule_id']]['source']=record['effective_source']
    for record in bundle['rules']:
        rule=record['template']
        if record['origin']!='book_appendix':continue
        rid=record['rule_id'];delta=rule['transform'].get('delta_formula')
        neutral=isinstance(delta,dict) and not rule['transform'].get('intrinsic_charge')
        source=bundle['pack_metadata'][record['compatibility_pack']]['source']
        by_id[rid]=dict(rule_id=rid,name_cn=rule['product_description'],origin='book_appendix',
            source={**source,**rule['source_refs']},mass_reference=rule['reference_state'],
            delta_formula=delta if neutral else None,delta_mass_mono_da=formula_mass(delta) if neutral else None,
            mass_status='formula_confirmed' if neutral and delta else 'isomeric' if neutral else 'structure_required',
            nominal_mass_shift_book=f"{rule['book_nominal_delta']:+} amu",
            diagnostic_evidence=rule['diagnostic_evidence'],relative_profile={
                'stage_tags':rule['required_conditions'],'any_residues':rule['site_selector'].get('residues',[]),
                'features':rule.get('required_features',[])})
        for alias in rule.get('legacy_rule_ids',[]):
            if alias in by_id:by_id[alias].update(alias_of=rid,automatic_mass_event=False,query_kind='legacy_alias')
    for record in bundle['rules']:
        facts=[x['content'] for x in record['review_layers'] if x['kind']=='facts']
        disposition=record['disposition']
        if disposition and disposition['status']=='compiled_fact':facts.insert(0,disposition.get('fact',{}))
        for fact in facts:
            by_id[record['rule_id']].update(query_kind='unresolved_phenomenon' if fact.get('reference_mass_lead') else 'condition_fact',
                automatic_mass_event=False,nominal_reference_lookup_allowed=fact.get('reference_mass_lead',False))
    return by_id,version


def add(*formulas):
    result = {}
    for formula in formulas:
        for element, count in formula.items():
            if element not in MONO or type(count) is not int:
                raise ValueError(f"Unsupported element or count: {element}={count}")
            result[element] = result.get(element, 0) + count
    return {k: v for k, v in result.items() if v}


def structure(record):
    return linear_formula(normalize_peptide(record))


def linear_formula(record):
    residues = record["sequence"]
    nterm, cterm = record["nterm"], record["cterm"]
    mods = record["modifications"]
    formula = add(parse_formula("H2O"), *(parse_formula(RESIDUE[r]) for r in residues),
                  {"C": 2, "H": 2, "O": 1} if nterm == "acetyl" else {},
                  C_TERM_DELTAS[cterm],
                  *(SULFUR_BOND_DELTAS[b["kind"]] for b in record.get("bonds", [])),
                  *(group_formula(m["group"]) if "group" in m else m["delta_formula"] for m in mods))
    if any(n < 0 for n in formula.values()):
        raise ValueError("Negative absolute element count")
    occupied = {m["site"] for m in mods} | {i for b in record.get("bonds", []) for i in b["sites"]}
    if any(m.get("group") in {"NCA","NChloroformylAcidChloride"} for m in mods):
        occupied.add("C")
    return residues, formula, occupied


def coordinate(observation):
    # Preserve the established CLI adapter's missing-field exception contract.
    if observation.get('kind')=='mz' and 'charge' not in observation:
        raise KeyError('charge')
    model=parse_coordinate(observation)
    return model.value, observation['tolerance_unit'], observation['tolerance'], model.multiplier, model.offset


def events(target, residues, occupied, rules, context):
    result = []

    def emit(rid, site, product, consumed):
        rule = rules[rid]
        if rule["mass_status"] not in {"formula_confirmed", "isomeric"}:
            raise ValueError(f"Activated rule {rid} has unsupported mass status")
        support, missing, conflict = stage_fit(rid, context)
        condition = "conflict" if conflict else "compatible" if support and not missing else "unknown"
        result.append({"rule_id": rid, "site": site, "product": product,
                       "consumed": list(consumed), "condition": condition,
                       "supported_conditions": support, "missing_conditions": missing,
                       "conflicting_conditions": conflict,
                       "delta_formula": rule["delta_formula"] or {},
                       "origin": rule["origin"], "source": rule.get("source", rule.get("source_basis")),
                       "diagnostic_evidence": rule.get("diagnostic_evidence", [])})

    for i, residue in enumerate(residues, 1):
        if i in occupied:
            continue
        for rid, (required, product) in SITE_RULES.items():
            if residue == required:
                emit(rid, i, product, (i,))
        if residue == "D" and i < len(residues):
            emit("SUP-003", i, "isoAsp", (i,))
            if i + 1 not in occupied and residues[i] != "P":
                emit("SUP-029", i, "Aspartimide", (i, i + 1))
    if target["nterm"] == "free" and not ({"N", 1} & occupied) and residues[0] in {"E", "Q"}:
        emit("SUP-004" if residues[0] == "E" else "SUP-005", 1, "pGlu", ("N", 1))
    if target["cterm"] == "amide" and "C" not in occupied:
        emit("SUP-045", "C", "C端羧酸", ("C",))
    for residue, group in itertools.groupby(enumerate(residues, 1), key=lambda item: item[1]):
        positions = tuple(i for i, _ in group)
        if residue != "P" or len(positions) < 2 or set(positions) & occupied:
            continue
        emit("CH10-R005", f"{positions[0]}–{positions[-1]}", "连续Pro区段多掺一枚Pro", positions)
        first = positions[0] - 1
        result[-1]["sequence_after_this_event"] = residues[:first] + ["P"] + residues[first:]
    return result


def diagnose(case, library=None):
    rules, version = library or load_rules()
    check_keys(case, {"case_id", "reaction", "context", "intended_product", "starting_material",
                      "observations", "max_events"}, "case")
    if "case_id" in case and (not isinstance(case["case_id"], str) or not case["case_id"].strip()):
        raise ValueError("case_id must be a nonempty display label")
    reaction = case["reaction"]
    check_keys(reaction, {"description", "stage_tags", "stage_tags_absent", "acyl_donors", "carbonyl_donors"}, "reaction")
    if not isinstance(reaction["description"], str) or not reaction["description"].strip():
        raise ValueError("Provide the actual intended reaction")
    tags = string_set(reaction.get("stage_tags", []), "stage_tags")
    if not isinstance(case["observations"], list) or not case["observations"] or not all(
        isinstance(item, dict) for item in case["observations"]
    ):
        raise ValueError("Provide at least one observation")
    target = normalize_peptide(case["intended_product"]) if case.get("intended_product") is not None else None
    starting = normalize_peptide(case["starting_material"]) if case.get("starting_material") is not None else None
    context = read_context(case, tags, target["sequence"] if target is not None else None)
    add_starting_group_facts(starting,target,context)
    donors = reaction.get("acyl_donors", [])
    if not isinstance(donors, list):
        raise ValueError("reaction.acyl_donors requires an explicit peptide list")
    context[0]["acyl_donors"] = [normalize_peptide(d) for d in donors]
    context[0]["carbonyl_donors"] = normalize_carbonyl_donors(reaction.get("carbonyl_donors", []))
    if any(d["cterm"] != "acid" or any(m["site"] == "C" for m in d["modifications"])
           or sequence_rule_exclusions(d) for d in context[0]["acyl_donors"]):
        raise ValueError("An acyl donor requires a defined free carboxylic acid and supported local groups")
    kinds = {item["kind"] for item in case["observations"]}
    relative_lookup = None
    if kinds & RELATIVE_KINDS:
        if not kinds <= RELATIVE_KINDS:
            raise ValueError("Use a separate case for absolute mass/mz and relative mass differences")
        appendix_source = (BOOK / "appendix_I_recovered.json").read_bytes()
        appendix_ledger = (BOOK / "appendix_I_element_account.json").read_bytes()
        ledger = {r["id"]: r for r in json.loads(appendix_ledger)["rows"]}
        appendix = [(r, ledger[r["appendix_row_id"]]) for r in json.loads(appendix_source)["rows"]]
        relative_lookup = list(search_relative(case, rules, context, appendix))
        if target is None or kinds == {"nominal_delta"}:
            if "max_events" in case:
                raise ValueError("Relative source lookup does not enumerate combinations")
            if any(item["reference"] == "starting_material" for item in case["observations"]):
                raise ValueError("Starting material comparison requires both the starting material and intended product")
            return {"engine_version": 5, "library_hash": version, "case": case,
                    "appendix_hash": hashlib.sha256(appendix_source + appendix_ledger).hexdigest(),
                    "observations": relative_lookup,
                    "scope": {"rules_indexed": len(rules), "appendix_rows_indexed": len(appendix),
                              "mode": "relative_source_lookup"}}
    if target is None:
        raise ValueError("Absolute mass requires intended_product")
    residues, formula, occupied = linear_formula(target)
    base = formula_mass(formula)
    start_formula = linear_formula(starting)[1] if starting is not None else None
    start_mass = formula_mass(start_formula) if start_formula is not None else None
    depth = case.get("max_events", 1)
    if type(depth) is not int or depth not in (1, 2):
        raise ValueError("max_events must be 1 or 2")
    singles = events(target, residues, occupied, rules, context)
    if (starting is not None and starting["sequence"] == target["sequence"]
            and all(starting[term] == target[term] for term in ("nterm", "cterm"))):
        for mod in starting.get("modifications", []):
            group, site = mod.get("group"), mod["site"]
            footprint = ["N", 1] if site == "N" else ["C", len(residues)] if site == "C" else [site]
            if group not in RETENTION_RULES or set(footprint) & occupied:
                continue
            rid = RETENTION_RULES[group]
            delta = group_formula(group)
            if delta != rules[rid]["delta_formula"]:
                raise ValueError(f"Protection group formula disagrees with {rid}")
            # N-terminal modifications can also interact with events on residue 1.
            support, missing, conflict = stage_fit(rid, context)
            condition = "conflict" if conflict else "unknown" if missing else "compatible"
            singles.append({"rule_id": rid, "site": site, "product": f"{group}残留",
                            "consumed": footprint, "condition": condition, "delta_formula": delta,
                            "supported_conditions": support, "missing_conditions": missing,
                            "conflicting_conditions": conflict,
                            "origin": rules[rid]["origin"], "source": rules[rid]["source_basis"],
                            "diagnostic_evidence": []})
    # State-only baseline migrations are executed by state_diagnosis, not this historical compatibility search.
    linear_raw = [path.read_bytes() for path in sorted((BOOK / "rules").glob("*_linear.json"))
                  if '_state_baselines_' not in path.name]
    bond_raw = [path.read_bytes() for path in sorted((BOOK / "rules").glob("*_bond_rules.json"))]
    linear_packs = [read_rule_pack(raw) for raw in linear_raw]
    linear_ids = [rule["rule_id"] for pack in linear_packs for rule in pack["rules"]]
    bond_packs = [read_rule_pack(raw) for raw in bond_raw]
    linear_ids += [rule["rule_id"] for pack in bond_packs for rule in pack["rules"]]
    if len(linear_ids) != len(set(linear_ids)) or not set(linear_ids) <= rules.keys():
        raise ValueError("Sequence rule packs require unique library identifiers")
    linear_events = [event for pack in linear_packs for event in instantiate(
        pack, target, context, linear_formula, add)]
    linear_events += [event for pack in bond_packs for event in instantiate_bond_rules(
        pack, target, context, linear_formula, add)]
    # A reviewed structural event replaces the legacy mass-only duplicate at that site.
    reviewed_sites = {(e["rule_id"], json.dumps(e["site"], sort_keys=True)) for e in linear_events}
    superseded = {e["rule_id"] for e in linear_events if e.get("supersedes_legacy")}
    singles = [e for e in singles if e["rule_id"] not in superseded
               and (e["rule_id"], json.dumps(e["site"], sort_keys=True)) not in reviewed_sites]
    independent = [e for e in linear_events if e.get("combination_policy") == "independent_other_sites"]
    singles += independent
    variants = [("intended_product", formula, [])]
    if start_formula is not None:
        variants.append(("starting_material", start_formula, []))
    for event in linear_events:
        if event.get("combination_policy") == "independent_other_sites":
            continue
        variants.append(("side_product_hypothesis", add(formula, event["delta_formula"]), [event]))
    for size in range(1, depth + 1):
        for combo in itertools.combinations(singles, size):
            sites = [site for e in combo for site in e["consumed"]]
            if len(sites) != len(set(sites)):
                continue
            composition = add(formula, *(e["delta_formula"] for e in combo))
            if all(n >= 0 for n in composition.values()):
                variants.append(("side_product_hypothesis", composition, list(combo)))

    observations = []
    active_ids = {e["rule_id"] for e in singles + linear_events}
    for index, observation in enumerate(case["observations"]):
        if observation["kind"] == "nominal_delta":
            observations.append(relative_lookup[index])
            continue
        if relative_lookup is None:
            measured, unit, tolerance, z, offset = coordinate(observation)
        else:
            ref = start_mass if observation["reference"] == "starting_material" else base
            if ref is None:
                raise ValueError("Starting material comparison requires starting_material")
            measured, unit, tolerance, z, offset = observation["value"], "Da", observation["tolerance"], 1, -ref

        def compare(mass):
            predicted = (mass + offset) / z
            error = measured - predicted
            if relative_lookup is not None:
                return abs(error) <= tolerance, predicted, error, None
            ppm = 1e6 * error / predicted if predicted > 0 else None
            match = predicted > 0 and abs(error if unit == "Da" else ppm) <= tolerance
            return match, predicted, error, ppm

        hits, conflicts = [], []
        basis=observation.get('mass_basis','average' if observation['kind']=='neutral_average' else 'monoisotopic')
        for kind, composition, combo in variants:
            intrinsic = sum(e.get("intrinsic_charge", 0) for e in combo)
            if intrinsic != observation.get("intrinsic_charge", 0):
                continue
            mass = formula_mass(composition)
            match, predicted, error, ppm = compare(observation_formula_mass(composition,basis))
            if match:
                bucket = conflicts if any(e["condition"] == "conflict" for e in combo) else hits
                bucket.append({"kind": kind, "events": combo, "formula": composition,
                             "neutral_mono": mass if not intrinsic else None,
                             **({'neutral_average':observation_formula_mass(composition,'average'),
                                 'observation_mass_basis':'average'} if basis=='average' else {}),
                             "delta_vs_target": mass - base if not intrinsic else None, "predicted": predicted,
                             **({"intrinsic_charge": intrinsic, "intrinsic_ion_mono": mass - intrinsic * ELECTRON,
                                 "element_mass_delta_vs_target": mass - base} if intrinsic else {}),
                             "err_da": error, "err_ppm": ppm, "evidence": "mass_match_only"})
        hits.sort(key=lambda h: (bool(h["events"]), sum(e["condition"] == "unknown" for e in h["events"]),
                                 len(h["events"]), abs(h["err_da"])))
        if relative_lookup is not None:
            result = relative_lookup[index]
            computed = {event["rule_id"] for hit in hits + conflicts for event in hit["events"]}
            result["reference_only_delta_leads"] = [
                lead for lead in result["reference_only_delta_leads"] if lead["rule_id"] not in computed]
        else:
            leads = []
            for rid, rule in rules.items():
                if observation.get("intrinsic_charge") or rid in REVIEWED_IDS or rid in active_ids or not has_neutral_lookup_delta(rule):
                    continue
                # The rule's substrate and mass reference are unverified here.
                delta = rule.get("delta_mass_mono_da")
                lead_base=base
                if basis=='average':
                    lead_base=observation_formula_mass(formula,'average')
                    delta=observation_formula_mass(rule['delta_formula'],'average') if isinstance(rule.get('delta_formula'),dict) else 0 if rule.get('mass_status')=='isomeric' else None
                if isinstance(delta, (int, float)) and compare(lead_base + delta)[0]:
                    lead = {k: rule[k] for k in (
                        "rule_id", "name_cn", "origin", "source", "source_basis", "mass_reference",
                        "susceptible_sequences", "applies_to", "conditions", "competing_explanations",
                        "diagnostic_evidence", "delta_mass_mono_da") if k in rule}
                    fit, support, missing, conflicting = applicability(rid, context,rule.get("relative_profile"))
                    leads.append(dict(lead, context_status=fit, support=support, missing=missing, conflicts=conflicting))
            result = {"input": observation, "reference_only_delta_leads": leads, "status": "no_automatic_match"}
        result.update(matches=hits, condition_conflicts=conflicts)
        if hits:
            result["status"] = "matched_hypotheses"
        observations.append(result)
    return {"engine_version": 5, "library_hash": version, "case": case, "start_neutral_mono": start_mass,
            **({"appendix_hash": hashlib.sha256(appendix_source + appendix_ledger).hexdigest()}
               if relative_lookup is not None else {}),
            "sequence_rule_hash": hashlib.sha256(b"\n".join(linear_raw+bond_raw)).hexdigest(),
            "target_neutral_mono": base, "intended_delta": base - start_mass if start_mass is not None else None,
            "observations": observations,
            "scope": {"rules_indexed": len(rules), "reviewed_rules": len(REVIEWED_IDS),
                      "active_for_case": sorted(active_ids),
                      "sequence_rule_templates": len(linear_ids),
                      "sequence_rules_skipped": bool(sequence_rule_exclusions(target)),
                      "sequence_rules_skip_reason": ("需要连接/构型或未定修饰重映射的普通模板未执行；"
                                                     "已明确的专用模板仍按声明范围运行，不能据此排除其他副产物。"
                                                     if sequence_rule_exclusions(target) else None),
                      "sequence_rules_blocking_modifications": sequence_rule_exclusions(target),
                      "component_dependent_rules_skipped": [r["rule_id"] for p in linear_packs for r in p["rules"]
                          if r["transform"]["kind"] in {"o_acylate", "acyl_donor_self_condensation"}
                          and not context[0]["acyl_donors"]],
                      "carbonyl_component_rules_skipped": [r["rule_id"] for p in linear_packs for r in p["rules"]
                          if r["transform"]["kind"]=="carbonyl_condensation" and not context[0]["carbonyl_donors"]],
                      "combinations": "independent legacy and explicitly declared local sites; topology/sequence edits single-event",
                      "max_events": depth}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = diagnose(json.loads(args.case.read_text(encoding="utf-8-sig")))
    except (ValueError, KeyError, TypeError) as error:
        parser.error(str(error))
    output = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)
    if args.output:
        with args.output.open("x", encoding="utf-8") as file:
            file.write(output + "\n")
    else:
        print(output)


if __name__ == "__main__":
    main()
