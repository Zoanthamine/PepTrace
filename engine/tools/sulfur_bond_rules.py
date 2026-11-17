"""Apply source-reviewed transformations to explicitly identified sulfur bridges."""
from copy import deepcopy
import itertools
import math
from case_input import evaluate_requirements


def disulfide_networks(target, limit):
    """Enumerate all matchings with the declared SS count; never silently truncate."""
    ss = [b for b in target.get("bonds", []) if b["kind"] == "disulfide"]
    count = len(ss)
    if not count:
        return
    fixed = [b for b in target["bonds"] if b["kind"] != "disulfide"]
    unavailable = {m["site"] for m in target["modifications"]} | {i for b in fixed for i in b["sites"]}
    sites = [i for i, r in enumerate(target["sequence"], 1) if r == "C" and i not in unavailable]
    possible = math.comb(len(sites), count*2) * math.factorial(count*2) // (2**count * math.factorial(count))
    if possible > limit:
        raise ValueError(f"Disulfide network enumeration requires {possible} states; declared rule limit is {limit}; no result was silently truncated")
    original = {tuple(b["sites"]) for b in ss}

    def perfect(vertices):
        if not vertices:
            yield []
            return
        first = vertices[0]
        for j, second in enumerate(vertices[1:], 1):
            for rest in perfect(vertices[1:j]+vertices[j+1:]):
                yield [(first, second)]+rest

    for chosen in itertools.combinations(sites, count*2):
        for pairs in perfect(list(chosen)):
            if set(pairs) == original:
                continue
            product = deepcopy(target)
            product["bonds"] = deepcopy(fixed)+[{"kind": "disulfide", "sites": list(p)} for p in pairs]
            yield product


def instantiate_bond_rules(pack, target, context, formula_of, add):
    parent = formula_of(target)[1]
    present, absent = context
    for rule in ({**base, **branch} for base in pack["rules"]
                 for branch in [{}]+base.get("alternate_branches", [])):
        if rule["transform"]["kind"] == "scramble_disulfides":
            support, missing, conflicts = evaluate_requirements(rule["required_conditions"], rule.get("required_any_conditions", []),
                                                                present["stage_tags"], absent["stage_tags"])
            for product in disulfide_networks(target, rule["transform"]["max_network_states"]):
                delta = add(formula_of(product)[1], {e: -n for e, n in parent.items()})
                if delta:
                    raise ValueError("Equal-count SS network changed elemental formula")
                yield dict(rule_id=rule["rule_id"], site=[b["sites"] for b in product["bonds"]],
                           consumed=[i for b in product["bonds"] for i in b["sites"]], delta_formula={},
                           product=rule["product_description"], product_structure=product,
                           event_kind="network_state_hypothesis", elementary_exchange_steps="not_enumerated",
                           condition="conflict" if conflicts else "unknown" if missing else "compatible",
                           supported_conditions=support, missing_conditions=missing, conflicting_conditions=conflicts,
                           origin="book", source={**pack["source"], **rule["source_refs"]},
                           source_evidence_level=rule["evidence_level"], reference_state=rule["reference_state"],
                           diagnostic_evidence=rule["diagnostic_evidence"])
            continue
        for index, bond in enumerate(target.get("bonds", [])):
            if bond["kind"] not in rule["bond_selector"]["kinds"]:
                continue
            first, second = bond["sites"]
            transform = rule["transform"]
            kind = transform["kind"]
            if kind == "oxidize_sulfur_bridge":
                positions = [first, second] if bond["kind"] == "disulfide" else [bond["oxidized_site"]]
            elif kind in {"disulfide_beta_elimination", "disulfide_alpha_elimination", "hydrolyze_disulfide"}:
                positions = [first, second]
            elif kind == "convert_sulfur_bridge" and transform["product_kind"] == "lanthionine":
                positions = [first, second]
            else:
                positions = [None]
            for site in positions:
                product = deepcopy(target)
                product_bond = product["bonds"][index]
                if kind == "oxidize_sulfur_bridge":
                    product_bond.update(kind=transform["product_kind"], oxidized_site=site)
                elif kind == "reduce_sulfur_bridge":
                    product["bonds"].pop(index)
                elif kind == "convert_sulfur_bridge":
                    product_bond["kind"] = transform["product_kind"]
                    product_bond.pop("oxidized_site", None)
                    if transform["product_kind"] == "lanthionine":
                        product.get("alpha_stereochemistry", {}).pop(str(site), None)
                elif kind in {"hydrolyze_disulfide", "disulfide_alpha_elimination"}:
                    product["bonds"].pop(index)
                    product["modifications"].append({"site": site, "group": "CysSulfenic" if kind == "hydrolyze_disulfide" else "CysThioaldehyde"})
                elif kind == "oxidative_bridge_cleavage":
                    product["bonds"].pop(index)
                    product["modifications"] += [{"site": i, "group": "CysSulfonic"} for i in bond["sites"]]
                elif kind == "hydrolyze_thiosulfinate":
                    product["bonds"].pop(index)
                    product["modifications"] += [{"site": i, "group": "CysSulfenic"} for i in bond["sites"]]
                elif kind == "disulfide_beta_elimination":
                    product["bonds"].pop(index)
                    partner = second if site == first else first
                    product["sequence"][site-1] = "Dha"
                    product.get("alpha_stereochemistry", {}).pop(str(site), None)
                    product["modifications"].append({"site": partner, "group": "Persulfide"})
                else:
                    raise ValueError(f"Unsupported bond transform: {kind}")
                if not product["bonds"]:
                    product.pop("bonds")
                    product.pop("topology", None)
                if "alpha_stereochemistry" in product and not product["alpha_stereochemistry"]:
                    product.pop("alpha_stereochemistry")
                delta = add(formula_of(product)[1], {e: -n for e, n in parent.items()})
                if "delta_formula" in transform and delta != transform["delta_formula"]:
                    raise ValueError(f'Inconsistent bond delta for {rule["rule_id"]}')
                support, missing, conflicts = evaluate_requirements(rule["required_conditions"],
                    rule.get("required_any_conditions", []), present["stage_tags"], absent["stage_tags"])
                yield dict(rule_id=rule["rule_id"], site=bond["sites"], consumed=bond["sites"], delta_formula=delta,
                           product=rule["product_description"], product_structure=product,
                           bond_edit={"input": bond, "reacting_site": site, "operation": kind},
                           condition="conflict" if conflicts else "unknown" if missing else "compatible",
                           supported_conditions=support, missing_conditions=missing, conflicting_conditions=conflicts,
                           origin="book", source={**pack["source"], **rule["source_refs"]},
                           source_evidence_level=rule["evidence_level"], reference_state=rule["reference_state"],
                           diagnostic_evidence=rule["diagnostic_evidence"])
