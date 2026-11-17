"""Instantiate source-reviewed sequence rules on explicitly linear peptides."""
import json

from case_input import evaluate_requirements
from mass_catalog import group_formula


def read_rule_pack(raw):
    pack = json.loads(raw)
    if pack["schema_version"] != 1:
        raise ValueError("Unsupported sequence rule schema")
    ids = [rule["rule_id"] for rule in pack["rules"]]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate sequence rule identifiers")
    return pack


def selected_sites(selector, residues):
    motif = selector.get("motif", ())
    position = selector["position"]
    if selector.get("sequence_length") and len(residues)!=selector["sequence_length"]:
        return
    if not set(selector.get("requires_residues", [])) <= set(residues):
        return
    for index, residue in enumerate(residues):
        if selector.get("next_residues") and (
            index + 1 == len(residues) or residues[index+1] not in selector["next_residues"]
        ):
            continue
        if selector.get("next_is_end") and index + 2 != len(residues):
            continue
        if selector.get("next_backbone_amide") and (
            index + 1 == len(residues) or residues[index + 1] == "P"
        ):
            continue
        if selector.get("previous_required") and index == 0:
            continue
        if position == "start" and index != 0:
            continue
        if position == "end" and index != len(residues) - 1:
            continue
        if position == "internal" and index in (0, len(residues) - 1):
            continue
        if position == "internal_or_end" and index == 0:
            continue
        if "residues" in selector and residue not in selector["residues"]:
            continue
        if residue in selector.get("exclude_residues", ()):
            continue
        if motif and residues[index:index + len(motif)] != motif:
            continue
        yield index


def slice_peptide(peptide, start, stop):
    end = len(peptide["sequence"])
    mods = []
    for mod in peptide["modifications"]:
        site = mod["site"]
        if type(site) is int and start < site <= stop:
            mods.append(dict(mod, site=site - start))
        elif (site == "N" and start == 0) or (site == "C" and stop == end):
            mods.append(mod)
    return {"sequence": peptide["sequence"][start:stop],
            "nterm": peptide["nterm"] if start == 0 else "free",
            "cterm": peptide["cterm"] if stop == end else "acid",
            "modifications": mods,
            **({'_origin_sites':peptide['_origin_sites'][start:stop]} if '_origin_sites' in peptide else {})}

def inserted_origin_sites(target, rule, cut, count):
    if '_origin_sites' not in target:return {}
    origins=target['_origin_sites'];anchor=origins[cut] if cut<len(origins) else origins[-1]+':C'
    prefix=anchor+'@'+rule.get('_branch_id',rule['rule_id'])
    existing=sum(x.startswith(prefix+'#') for x in origins)
    added=[prefix+f'#{existing+i+1}' for i in range(count)]
    return {'_origin_sites':origins[:cut]+added+origins[cut:]}


def linked_products(target, idx, sites, transform):
    """Generate same-peptide linkage candidates and retain chain multiplicity."""
    bond = transform["bond"]
    symmetric = transform.get("symmetric", False)
    for partner in sites:
        if not symmetric or partner >= idx:
            yield ("同序列双链", {
                "topology": "crosslinked", "chains": [target, target],
                "bond": {"kind": bond, "sites": [[1, idx + 1], [2, partner + 1]]},
            }, transform["bond_delta"])
        if transform.get("intrachain") and idx < partner:
            yield ("单链内交联", {
                "topology": "crosslinked", "chains": [target],
                "bond": {"kind": bond, "sites": [[1, idx + 1], [1, partner + 1]]},
            }, transform["bond_delta"])


def sequence_rule_exclusions(target):
    """Explain when sequence edits cannot safely preserve a modification."""
    unsupported = {"Disulfide", "HeadTailCyclization", "Hydrolysis", "WaterLoss", "AmmoniaLoss", "NCA", "NChloroformylAcidChloride"}
    blocked = [dict(mod) for mod in target["modifications"]
            if mod.get("group") in unsupported or ("delta_formula" in mod and not
               (target.get('_atom_model') and target.get('_origin_sites') and ('atom' in mod or 'atoms' in mod)))] + [
                {"bond": dict(b), "reason": "backbone edits on sulfur-linked inputs require topology preservation"}
                for b in target.get("bonds", [])]
    if target.get("alpha_stereochemistry") and not (target.get('_atom_model') and target.get('_origin_sites')):
        blocked.append({"stereochemistry": target["alpha_stereochemistry"],
                        "reason": "general sequence edits need stereocenter remapping"})
    return blocked


def bridge_products(target, idx, transform, occupied):
    """Enumerate explicit atom-labelled bridges on one chain or identical chains.

    This describes connectivity and stoichiometry; it does not assert a conformation.
    Separate partner atoms and external free molecules stay separate candidates.
    """
    first = transform.get("source_site", idx + 1)
    allowed_groups = set(transform.get("allowed_occupied_groups", []))
    first_allowed = any(m["site"] == first and m.get("group") in allowed_groups for m in target["modifications"])
    if (first in occupied and not first_allowed) or (first == "N" and target["nterm"] != "free"):
        return
    external = transform.get("external_partner")
    partner_chain = external or target
    partner_occupied = set() if external else occupied
    for atom, selector in transform["partner_atoms"].items():
        for j in selected_sites(selector, partner_chain["sequence"]):
            site = selector.get("site_label", j + 1)
            partner_allowed = any(m["site"] == site and m.get("group") in allowed_groups for m in partner_chain["modifications"])
            required_group = selector.get("group")
            if required_group and not any(m["site"] == site and m.get("group") == required_group for m in partner_chain["modifications"]):
                continue
            if (site in partner_occupied and not partner_allowed) or (site == "N" and partner_chain["nterm"] != "free"):
                continue
            left = {"chain": 1, "site": first, "atom": transform["source_atom"]}
            for intra in ([False] if external else [True, False]):
                if transform.get("topologies") and ("intrachain" if intra else "interchain") not in transform["topologies"]:
                    continue
                if intra and (site == first or transform.get("symmetric") and j <= idx):
                    continue
                if not intra and transform.get("symmetric") and j < idx:
                    continue
                right = {"chain": 1 if intra else 2, "site": site, "atom": atom}
                chains = [target] if intra else [target, partner_chain]
                yield ("链内桥" if intra else "外加游离分子桥" if external else "同序列双链桥",
                       {"topology": "crosslinked", "chains": chains,
                        "bonds": [{"kind": transform["bond"], "endpoints": [left, right],
                                   "delta_formula": transform["bond_delta"],
                                   **({"consumed_unsaturation": {"chain": 1, "site": first, "kind": transform["source_unsaturation"]},
                                       "new_alpha_configuration": "not_assigned"} if transform.get("source_unsaturation") else {})}],
                        "geometry_status": "not_assessed"}, transform["bond_delta"])


def instantiate(pack, target, context, formula_of, add):
    """Yield each product independently, with its own formula and site identity."""
    residues, parent, occupied = formula_of(target)
    exclusions = sequence_rule_exclusions(target)
    present, absent = context
    branches = ({**rule, **branch}
                for rule in pack["rules"] for branch in [{}] + rule.get("alternate_branches", []))
    for rule in branches:
        selector, transform = rule["site_selector"], rule["transform"]
        composite_conversion = (transform["kind"]=="replace_group" and selector.get("group")=="NCA"
                                and all(item.get("group")=="NCA" for item in exclusions))
        if exclusions and not composite_conversion and (transform["kind"] != "epimerize_alpha"
                           or any("bond" not in item and "stereochemistry" not in item for item in exclusions)):
            continue
        if "cterm" in selector and selector["cterm"] != target["cterm"]:
            continue
        if "nterm" in selector and selector["nterm"] != target["nterm"]:
            continue
        eligible = list(selected_sites(selector, residues))
        for idx in eligible:
            kind = transform["kind"]
            if selector.get("nterm_participates") and (target["nterm"]!="free" or "N" in occupied):
                continue
            count = transform.get("count", 1)
            positions = list(range(idx + 1, idx + 1 + count))
            if kind == "duplicate_residue":
                positions = [idx + 1]
            if selector.get("next_backbone_amide"):
                positions.append(idx + 2)
            if kind == "c_terminal_amide_substitution":
                positions = ["C"]
            if kind == "n_terminal_delta" and selector.get("side_chain_required"):
                if idx + 1 in occupied:
                    continue
            group = selector.get("group")
            group_site = selector.get("group_site", idx + 1)
            if selector.get("groups"):
                group = next((mod.get("group") for mod in target["modifications"]
                              if mod["site"] == group_site and mod.get("group") in selector["groups"]), None)
                if group is None:
                    continue
            if selector.get("side_chain_required") and group_site == "N" and idx + 1 in occupied:
                continue
            if group is not None and not any(
                mod["site"] == group_site and mod.get("group") == group for mod in target["modifications"]
            ):
                continue
            bound_mods=[mod for mod in target['modifications'] if mod['site']==group_site and mod.get('group')==group] if group else []
            if group and target.get('_atom_model'):
                from atom_model import modification_matches_selector
                bound_mods=[m for m in bound_mods if modification_matches_selector(rule,target,m,group_site)]
                if not bound_mods:continue
            if group is None and kind not in {"n_terminal_delta", "capped_suffix", "site_bridge", "epimerize_alpha"} and occupied.intersection(positions):
                continue
            products = []
            stereo_edit = None
            if kind == "epimerize_alpha":
                if selector.get("unmodified_side_chain") and idx+1 in occupied:
                    continue
                product = dict(target)
                if transform.get("replace_residue"):
                    product["sequence"] = list(residues)
                    product["sequence"][idx] = transform["replace_residue"]
                before = target.get("alpha_stereochemistry", {}).get(str(idx+1))
                invert = transform.get("invert", True)
                after = {"L": "D", "D": "L"}.get(before) if invert else before
                stereo_edit = {"site": idx+1, "center": "alpha", "operation": "invert_relative_to_input" if invert else "retain_relative_to_input",
                               "input_LD": before, "product_LD": after,
                               "absolute_assignment": "declared_relative_LD" if before else "not_assigned"}
                if before:
                    product["alpha_stereochemistry"] = dict(target["alpha_stereochemistry"], **{str(idx+1): after})
                if transform.get("linkage"):
                    product["linkage_edit"] = {"site": idx+1, "kind": transform["linkage"]}
                products = [("相对原输入Cα构型反转" if invert else "保留原输入Cα构型", product, {})]
            elif kind == "hydrolytic_split":
                cut = idx if transform.get("before_site") else idx + 1
                if not 0 < cut < len(residues) or {cut, cut + 1} & occupied:
                    continue
                positions = [cut, cut + 1]
                products = [("N端片段", slice_peptide(target, 0, cut), {}),
                            ("C端片段", slice_peptide(target, cut, len(residues)), {})]
            elif kind == "backbone_amide_to_ester":
                if not idx or idx in occupied:
                    continue
                positions = [idx, idx + 1]
                products = [("同式酯异构体", dict(target, bond_edit={
                    "remove": [[idx, "carbonyl_C"], [idx + 1, "backbone_N"]],
                    "add": [[idx, "carbonyl_C"], [idx + 1, "side_chain_O"]],
                }), {})]
            elif kind == "n_acetyl_change":
                mods = target["modifications"]
                if transform["move_to_oxygen"]:
                    mods = mods + [{"site": idx + 1, "group": "Ac", "name": "侧链O-乙酰"}]
                products = [("乙酰转移产物" if transform["move_to_oxygen"] else "游离N端肽",
                             dict(target, nterm="free", modifications=mods), {})]
            elif kind == "n_terminal_delta":
                if "N" in occupied:
                    continue
                positions = ["N", idx + 1] if selector.get("side_chain_required") else ["N"]
                if selector.get("next_backbone_amide"):
                    positions.extend([idx+1,idx+2])
                if selector.get("terminal_amide_participates"):
                    positions.extend([idx+1,"C"])
                mods = target["modifications"] + [{"site": "N", **({"group":transform["product_group"]} if transform.get("product_group")
                                                                  else {"delta_formula":transform["delta_formula"]}),
                                                   "name": transform["modification"]}]
                products = [("N端修饰肽", dict(target, modifications=mods), {})]
            elif kind == "carbonyl_condensation":
                if target["nterm"]!="free" or "N" in occupied:
                    continue
                positions=["N"]
                for component in present.get("carbonyl_donors",[]):
                    delta=add(component["formula"],{"H":-2,"O":-1})
                    for carbonyl in component["carbonyl_sites"]:
                        configurations=["E","Z"] if carbonyl["distinct_substituents"] else ["not_applicable"]
                        for configuration in configurations:
                            mods=target["modifications"]+[{"site":"N","delta_formula":delta,"name":f"Oxime:{component['name']}:{configuration}"}]
                            product=dict(target,modifications=mods,condensation={
                                "kind":"oxime","component":component,"carbonyl_atom_index":carbonyl["atom_index"],
                                "oxime_configuration":configuration})
                            products.append((f"{component['name']}羰基{carbonyl['atom_index']}：{configuration}",product,{}))
            elif kind == "capped_suffix":
                # A synthesis interrupted at index idx retains its C-terminal fragment.
                # idx is in the full target's coordinates; the cap is on the new terminus.
                if idx == 0 and (target["nterm"] != "free" or "N" in occupied):
                    continue
                fragment = slice_peptide(target, idx, len(residues))
                mod = {"site": "N", "delta_formula": transform["cap_formula"],
                       "name": transform["modification"]}
                positions = ["N"] if not idx else list(range(1, idx + 1)) + [idx + 1]
                branch = "封端全长肽" if not idx else "封端C端截短片段"
                products = [(branch, dict(fragment, modifications=fragment["modifications"] + [mod]), {})]
            elif kind == "prefix_dkp":
                if idx != 0 or len(residues) <= count or target["nterm"] != "free" or "N" in occupied:
                    continue
                products = [("C端短肽", slice_peptide(target, count, len(residues)), {}),
                            ("环二肽", slice_peptide(target, 0, count), {"H": -2, "O": -1})]
            elif kind == "delete":
                if (idx == 0 and not transform.get("allow_n_terminal")) or idx + count >= len(residues):
                    continue
                if idx == 0 and (target["nterm"] != "free" or "N" in occupied):
                    continue
                mods = [dict(mod, site=mod["site"] - count)
                        if type(mod["site"]) is int and mod["site"] > idx + count else mod
                        for mod in target["modifications"]
                        if not (type(mod["site"]) is int and idx < mod["site"] <= idx+count)]
                products = [("缺失肽", dict(target, sequence=residues[:idx] + residues[idx+count:],
                                          modifications=mods,
                                          **({'_origin_sites':target['_origin_sites'][:idx]+target['_origin_sites'][idx+count:]} if '_origin_sites' in target else {})), {})]
            elif kind == "duplicate_residue":
                inserted = [residues[idx]] * count
                positions = [idx + 1]
                mods = [dict(mod, site=mod["site"] + count)
                        if type(mod["site"]) is int and mod["site"] > idx + 1 else mod
                        for mod in target["modifications"]]
                if transform.get("insert_phosphorylated"):
                    mods = mods + [{"site": idx + 2 + j, "group": "Phospho"} for j in range(count)]
                products = [("重复掺入肽", dict(target,
                    sequence=residues[:idx+1] + inserted + residues[idx+1:],
                    modifications=mods,**inserted_origin_sites(target,rule,idx+1,count)), {})]
            elif kind == "insert_gly_before":
                for cut in range(idx + 1):
                    mods = [dict(mod, site=mod["site"] + 1)
                            if type(mod["site"]) is int and mod["site"] > cut else mod
                            for mod in target["modifications"]]
                    inserted = dict(target, sequence=residues[:cut] + ["G"] + residues[cut:],
                                    modifications=mods,**inserted_origin_sites(target,rule,cut,1))
                    products.append((f"原序列位置{cut}与{cut+1}之间插入Gly", inserted, {}))
            elif kind == "insert_before":
                mods = [dict(mod, site=mod["site"] + 1)
                        if type(mod["site"]) is int and mod["site"] > idx else mod
                        for mod in target["modifications"]]
                products = [("掺入非天然残基", dict(target,
                    sequence=residues[:idx] + [transform["residue"]] + residues[idx:],
                    modifications=mods,**inserted_origin_sites(target,rule,idx,1)), {})]
            elif kind == "substitute":
                replacement = (transform["residues_by_input"][residues[idx]]
                               if "residues_by_input" in transform else transform["residue"])
                mods = [mod for mod in target["modifications"] if mod not in bound_mods or group is None]
                products = [("替换肽", dict(target, sequence=residues[:idx] + [replacement] +
                                             residues[idx+1:], modifications=mods), {})]
            elif kind == "local_delta":
                modification = {"site": idx + 1, **({"group": transform["product_group"]}
                    if transform.get("product_group") else {"delta_formula": transform["delta_formula"]}),
                                "name": transform["modification"]}
                products = [("局部修饰肽", dict(target, modifications=target["modifications"] + [modification]), {})]
            elif kind == "oxidize_protected_thioether":
                if group is None:
                    continue
                mods = [m for m in target["modifications"] if m not in bound_mods]
                mods.append({"site": group_site, "delta_formula": add(group_formula(group), {"O": 1}),
                             "name": f"Cys({group}) sulfoxide"})
                products = [("保护基保留的Cys sulfoxide", dict(target, modifications=mods), {})]
            elif kind == "remove_group":
                positions = [group_site]
                mods = [mod for mod in target["modifications"] if mod not in bound_mods]
                products = [("脱除基团", dict(target, modifications=mods), {})]
            elif kind == "transfer_group_to_amine":
                # Move a covalent group; elemental content is preserved.
                # Keep every distinct donor/acceptor interpretation.
                mods = [mod for mod in target["modifications"] if mod not in bound_mods]
                destinations = (["N"] if target["nterm"] == "free" and "N" not in occupied else [])
                destinations += [j + 1 for j, r in enumerate(residues)
                                 if r == "K" and j + 1 not in occupied]
                if transform.get("destination") == "N":
                    destinations = [s for s in destinations if s == "N"]
                for dest in destinations:
                    if dest == group_site:
                        continue
                    products.append((f"{group_site}→{dest}基团转移", dict(target, modifications=mods + [
                        {"site": dest, "group": group, "name": transform["modification"]}]), {}))
                positions = [group_site]
            elif kind == "transfer_group_to_oxygen":
                if idx + 1 in occupied:
                    continue
                mods = [mod for mod in target["modifications"] if mod not in bound_mods]
                positions = [group_site, idx + 1]
                products = [("N→O酰基转移", dict(target, modifications=mods + [
                    {"site": idx + 1, "group": group, "name": transform["modification"]}]), {})]
            elif kind == "c_terminal_delta":
                if "C" in occupied:
                    continue
                positions = ["C"]
                products = [("C端官能团修饰", dict(target, modifications=target["modifications"] + [
                    {"site": "C", "delta_formula": transform["delta_formula"],
                     "name": transform["modification"]}]), {})]
            elif kind == "replace_group":
                positions = [group_site, idx + 1] if selector.get("side_chain_required") else [group_site]
                if selector.get("nterm_participates"):
                    positions=["N",group_site]
                mods = [mod for mod in target["modifications"] if mod not in bound_mods]
                mods.append({"site": group_site, **({"group": transform["product_group"]} if transform.get("product_group")
                                                    else {"delta_formula": transform["replacement_formula"]}),
                             "name": transform["modification"]})
                products = [("转化保护基或修饰", dict(target, modifications=mods), {})]
            elif kind == "paired_group_redistribution":
                partners = [m["site"] for m in target["modifications"] if m.get("group") == group and type(m["site"]) is int and m["site"] > idx+1]
                for partner in partners:
                    mods = [m for m in target["modifications"] if m["site"] not in {idx+1, partner}]
                    for changed in [idx+1, partner]:
                        new = mods + [{"site": changed, "group": transform["product_group"]}]
                        products.append((f"两个{group}歧化：{changed}为{transform['product_group']}",
                                         dict(target, modifications=new, _event_sites=[idx+1, partner]), {}))
                positions = [idx+1]
            elif kind == "migrate_group_to_previous":
                if idx == 0 or idx in occupied or residues[idx - 1] != transform["previous_residue"]:
                    continue
                positions.insert(0, idx)
                moved = {"site": idx, "delta_formula": transform["migrated_formula"],
                         "name": transform["modification"]}
                mods = [mod for mod in target["modifications"] if mod not in bound_mods]
                products = [("基团转移", dict(target, modifications=mods + [moved]), {})]
            elif kind == "c_terminal_amide_substitution":
                modified = dict(target, cterm="n_4_hydroxybenzyl_amide")
                products = [("取代酰胺", modified, {})]
            elif kind == "c_terminal_pyrrolidinamide":
                if "C" in occupied:
                    continue
                positions = ["C"]
                products = [("吡咯烷酰胺", dict(target, cterm="pyrrolidinamide"), {})]
            elif kind == 'c_terminal_acid':
                if 'C' in occupied:continue
                positions=['C']
                products=[('C端酰胺水解为羧酸',dict(target,cterm='acid'),{})]
            elif kind == "trp_crosslink":
                products = linked_products(target, idx, [i for i in eligible if i + 1 not in occupied], transform)
            elif kind == "site_bridge":
                products = bridge_products(target, idx, transform, occupied)
            elif kind in {"o_acylate", "n_acylate"}:
                for donor in present.get("acyl_donors", []):
                    products.append(("显式酰基供体接至O/N", {"topology": "branched_acyl_peptide", "chains": [target, donor],
                        "bonds": [{"kind": "O_acyl_ester" if kind=="o_acylate" else "amide", "endpoints": [
                            {"chain": 1, "site": transform.get('source_site',idx + 1), "atom": "sidechain_O" if kind=="o_acylate" else transform['source_atom']},
                            {"chain": 2, "site": "C", "atom": "carbonyl_C"}],
                            "delta_formula": {"H": -2, "O": -1}}]}, {"H": -2, "O": -1}))
            elif kind == "acyl_donor_self_condensation":
                for donor in present.get("acyl_donors", []):
                    if not any(m["site"] == "N" and m.get("group") == "Fmoc" for m in donor["modifications"]):
                        continue
                    size = len(donor["sequence"])
                    if idx != size or residues[:size] != donor["sequence"]:
                        continue
                    mods = donor["modifications"] + [dict(m, site=m["site"] + size)
                        for m in donor["modifications"] if type(m["site"]) is int]
                    products.append(("供体去Fmoc后与另一供体自缩合，受体未入产物", dict(donor,
                        sequence=donor["sequence"]*2, modifications=mods), {}))
            else:
                raise ValueError(f"Unsupported sequence transform: {kind}")
            support, missing, conflicts = evaluate_requirements(rule["required_conditions"],
                rule.get("required_any_conditions", ()),
                present["stage_tags"], absent["stage_tags"])
            features = list(rule.get("required_features", ()))
            site_fact = selector.get("required_site_fact", rule.get("required_site_fact"))
            if site_fact:
                fact_site = selector.get("fact_site", idx + 1)
                features.append(f'{site_fact}:{fact_site}')
            if kind == "capped_suffix" and idx:
                features.append(f'{transform.get("boundary_fact", "synthesis_boundary")}:{idx+1}')
            yes, unknown, no = evaluate_requirements(features, (), present["features"], absent["features"])
            support.extend(yes)
            missing.extend(unknown)
            conflicts.extend(no)
            support.extend(sorted(set(rule.get("supporting_any_conditions", ())) & present["stage_tags"]))
            condition = "conflict" if conflicts else "compatible" if support and not missing else "unknown"
            for branch, product, topology_delta in products:
                if kind in {"trp_crosslink", "site_bridge", "o_acylate", "n_acylate"}:
                    composition = add(*(formula_of(chain)[1] for chain in product["chains"]), topology_delta)
                    if kind == "trp_crosslink":
                        event_site = product["bond"]["sites"]
                    else:
                        endpoints = product["bonds"][0]["endpoints"]
                        event_site = [[a["chain"], a["site"], a["atom"]] for a in endpoints]
                        positions = event_site
                    product_structure = product
                else:
                    if kind == "paired_group_redistribution":
                        positions = product.pop("_event_sites")
                    composition = add(formula_of(product)[1], topology_delta)
                    event_site = positions[0] if len(positions) == 1 else positions
                    topology = "depsipeptide" if kind == "backbone_amide_to_ester" else (
                        "head_tail_cycle" if topology_delta else "linear")
                    product_structure = dict(product, topology=topology)
                    if transform.get("local_cycle"):
                        product_structure.update(topology="local_cycle", cycle={
                            "kind": transform["local_cycle"], "site": idx + 1})
                    if transform.get("monomer_state"):
                        product_structure["monomer_state"]=transform["monomer_state"]
                    if kind == "transfer_group_to_amine":
                        dest = next(m["site"] for m in product["modifications"]
                                    if m.get("group") == group and m["site"] not in occupied)
                        event_site = [group_site, dest]
                        positions = [group_site, dest]
                delta = add(composition, {element: -n for element, n in parent.items()})
                if "delta_formula" in transform:
                    expected_delta=transform["delta_formula"]
                    if kind=='delete':
                        from knowledge_source import components
                        removed=[components()['covalent_groups'][m['group']]['net_formula'] if 'group' in m else m['delta_formula']
                            for m in target['modifications'] if type(m['site']) is int and idx<m['site']<=idx+count]
                        expected_delta=add(expected_delta,*[{e:-n for e,n in f.items()} for f in removed])
                    if delta != expected_delta:
                        raise ValueError(f'Rule {rule["rule_id"]} gives an inconsistent elemental change')
                yield {
                    "rule_id": rule["rule_id"], "site": event_site,
                    "product": f'{rule["product_description"]}：{branch}',
                    "product_structure": product_structure,
                    "condition": condition, "supported_conditions": support,
                    "missing_conditions": missing, "conflicting_conditions": conflicts,
                    "consumed": positions, "delta_formula": delta,
                    "origin": ("inference_from_book" if rule["evidence_level"] == "editorial_extrapolation_from_book"
                               else pack["source"].get("origin", "book")),
                    "source": {**pack["source"], **rule["source_refs"]},
                    "diagnostic_evidence": rule["diagnostic_evidence"], "source_evidence_level": rule["evidence_level"],
                    "reference_state": rule["reference_state"],
                    **({"stereochemistry_edit": stereo_edit} if stereo_edit else {}),
                    "combination_policy": rule.get("combination_policy", "single_event"),
                    "supersedes_legacy": rule.get("supersedes_legacy", False),
                    **({"intrinsic_charge": transform["intrinsic_charge"]} if "intrinsic_charge" in transform else {}),
                }
