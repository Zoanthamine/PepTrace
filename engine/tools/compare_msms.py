"""Compare linear b/y fragment hypotheses against one centroided MS2 spectrum."""
import argparse
from bisect import bisect_left, bisect_right
from importlib.metadata import version
import json
import math
from pathlib import Path
import re

from case_input import C_TERM_DELTAS, normalize_peptide
from mass_catalog import group_formula
from mass_ledger import formula_mass

ION = re.compile(r"([by])(\d+)(\++)$")
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
NONLOCAL = {"Disulfide", "HeadTailCyclization", "Hydrolysis", "WaterLoss", "AmmoniaLoss"}


def read_ms2(path, native_id):
    import pyopenms as oms

    run = oms.OnDiscMSExperiment()
    if not run.openFile(str(path)):
        raise ValueError("mzML needs a readable index; use a centroided peak list for unindexed files")
    scan = run.getSpectrumByNativeId(native_id)
    if scan.getNativeID() != native_id or scan.getMSLevel() != 2:
        raise ValueError("Requested native ID must identify one MS2 spectrum")
    if scan.getType() != oms.SpectrumSettings.SpectrumType.CENTROID:
        raise ValueError("Selected mzML spectrum must be labeled centroided")
    mzs, intensities = scan.get_peaks()
    return [[float(mz), float(intensity)] for mz, intensity in zip(mzs, intensities)], {
        "mzml": str(path), "native_id": native_id, "rt_seconds": scan.getRT(),
        "precursors": [{"mz": p.getMZ(), "charge": p.getCharge(),
                        "lower_offset": p.getIsolationWindowLowerOffset(),
                        "upper_offset": p.getIsolationWindowUpperOffset()} for p in scan.getPrecursors()],
    }


def fragments(peptide, max_charge):
    import pyopenms as oms

    peptide = normalize_peptide(peptide)
    residues = peptide["sequence"]
    if len(residues) < 2 or not set(residues) <= STANDARD:
        raise ValueError("Fragment model requires at least two standard residues")
    shifts = []
    for mod in peptide.get("modifications", []):
        if mod.get("group") in NONLOCAL:
            raise ValueError("A nonlocal change cannot be modeled as a residue shift")
        shift = group_formula(mod["group"]) if "group" in mod else mod["delta_formula"]
        shifts.append((mod["site"], formula_mass(shift)))
    if peptide["nterm"] == "acetyl":
        shifts.append(("N", formula_mass(group_formula("Ac"))))
    cterm_delta = C_TERM_DELTAS[peptide["cterm"]]
    if cterm_delta:
        shifts.append(("C", formula_mass(cterm_delta)))

    generator = oms.TheoreticalSpectrumGenerator()
    params = generator.getParameters()
    for key in ("add_metainfo", "add_first_prefix_ion", "add_b_ions", "add_y_ions"):
        params.setValue(key, "true")
    for key in ("add_losses", "add_term_losses", "add_precursor_peaks", "add_abundant_immonium_ions"):
        params.setValue(key, "false")
    params.setValue("isotope_model", "none")
    generator.setParameters(params)
    spectrum = oms.MSSpectrum()
    generator.getSpectrum(spectrum, oms.AASequence.fromString("".join(residues)), 1, max_charge)
    result = []
    for ion, mz in zip(spectrum.getStringDataArrays()[0], spectrum.get_peaks()[0]):
        series, length, charge = ION.fullmatch(ion).groups()
        length, charge = int(length), len(charge)
        if series == "b":
            sites = {"N", *range(1, length + 1)}
        else:
            sites = {"C", *range(len(residues) - length + 1, len(residues) + 1)}
        predicted = float(mz) + math.fsum(delta for site, delta in shifts if site in sites) / charge
        if predicted <= 0:
            raise ValueError(f"Nonpositive mass for fragment {ion}; check the assigned modifications")
        result.append((ion, predicted))
    return result


def compare(candidates, peaks, tolerance, unit, max_charge=1):
    if not isinstance(candidates, dict) or len(candidates) < 2:
        raise ValueError("Provide at least two named candidates")
    if type(max_charge) is not int or not 1 <= max_charge <= 6:
        raise ValueError("max_charge must be an integer in 1..6")
    if unit not in {"Da", "ppm"} or not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("Provide positive tolerance in Da or ppm")
    if not peaks:
        raise ValueError("Provide nonempty centroided MS2 peaks")
    indexed = []
    for i, (mz, intensity) in enumerate(peaks, 1):
        if not math.isfinite(mz) or not math.isfinite(intensity) or mz <= 0 or intensity < 0:
            raise ValueError("Peaks must be finite, with positive m/z and nonnegative intensity")
        indexed.append((mz, i, intensity))
    indexed.sort()
    mzs = [x[0] for x in indexed]
    annotations = {}
    for name, peptide in candidates.items():
        for ion, mz in fragments(peptide, max_charge):
            window = tolerance if unit == "Da" else tolerance * mz * 1e-6
            lo, hi = bisect_left(mzs, mz - window), bisect_right(mzs, mz + window)
            for observed, peak_id, intensity in indexed[lo:hi]:
                peak = annotations.setdefault(peak_id, {"peak_id": peak_id, "mz": observed,
                                                        "intensity": intensity, "matches": []})
                peak["matches"].append({"candidate": name, "ion": ion, "theoretical_mz": mz,
                                        "error_da": observed - mz, "error_ppm": (observed - mz) / mz * 1e6})
    matched = sorted(annotations.values(), key=lambda p: p["peak_id"])
    for peak in matched:
        names = {m["candidate"] for m in peak["matches"]}
        peak["exclusive_candidate_among_tested"] = next(iter(names)) if len(names) == 1 else None
    return {"model": "positive intact b/y ions with local elemental mass shifts", "backend": "pyopenms",
            "backend_version": version("pyopenms"), "candidates": candidates,
            "tolerance": tolerance, "tolerance_unit": unit, "max_fragment_charge": max_charge,
            "observed_peak_count": len(peaks), "matched_peaks": matched,
            "unmatched_peak_count": len(peaks) - len(matched)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--mzml", type=Path)
    parser.add_argument("--native-id")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        case = json.loads(args.case.read_text(encoding="utf-8-sig"))
        if bool(args.mzml) != bool(args.native_id) or (args.mzml and "peaks" in case):
            raise ValueError("Provide either case peaks or --mzml paired with --native-id")
        peaks, source = read_ms2(args.mzml, args.native_id) if args.mzml else (case["peaks"], {"peak_list": str(args.case)})
        result = compare(case["candidates"], peaks, case["tolerance"], case["tolerance_unit"],
                         case.get("max_fragment_charge", 1))
        result["spectrum_source"] = source
    except (ValueError, TypeError, KeyError, RuntimeError, ImportError, OSError) as error:
        parser.error(str(error))
    output = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)
    if args.output:
        with args.output.open("x", encoding="utf-8") as handle:
            handle.write(output + "\n")
    else:
        print(output)


if __name__ == "__main__":
    main()
