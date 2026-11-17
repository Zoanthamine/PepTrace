import json
from pathlib import Path
import tempfile
import unittest

import pyopenms as oms
from pyteomics import mass

from compare_msms import compare, fragments, read_ms2

EXAMPLE = json.loads(Path(__file__).with_name("example_msms_compare.json").read_text(encoding="utf-8"))


class FragmentTests(unittest.TestCase):
    def test_site_discrimination_and_common_ion(self):
        result = compare(EXAMPLE["candidates"], EXAMPLE["peaks"], 5, "ppm")
        self.assertEqual([1, 2, 3], [p["peak_id"] for p in result["matched_peaks"]])
        exclusive = [p["exclusive_candidate_among_tested"] for p in result["matched_peaks"]]
        self.assertEqual(["Met2+O", "Met2+O", None], exclusive)
        self.assertEqual(1, result["unmatched_peak_count"])

    def test_b_y_charge_and_termini_against_pyteomics(self):
        peptide = {"sequence": "GMMK", "nterm": "acetyl", "cterm": "amide",
                   "modifications": [{"site": 2, "group": "Oxidation"}]}
        ions = dict(fragments(peptide, 2))
        expected_b = mass.fast_mass("GM", ion_type="b", charge=2) + (42.010565 + 15.994915) / 2
        expected_y = mass.fast_mass("MK", ion_type="y", charge=2) - 0.984016 / 2
        self.assertAlmostEqual(ions["b2++"], expected_b, delta=0.00001)
        self.assertAlmostEqual(ions["y2++"], expected_y, delta=0.00001)

    def test_indexed_mzml_native_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "synthetic.mzML"
            experiment = oms.MSExperiment()
            scan = oms.MSSpectrum()
            scan.setNativeID("scan=25")
            scan.setMSLevel(2)
            scan.setType(oms.SpectrumSettings.SpectrumType.CENTROID)
            scan.set_peaks(([205.064139, 278.153289], [100, 200]))
            experiment.addSpectrum(scan)
            oms.MzMLFile().store(str(path), experiment)
            peaks, metadata = read_ms2(path, "scan=25")
            self.assertEqual("scan=25", metadata["native_id"])
            self.assertEqual(2, len(peaks))
            self.assertAlmostEqual(205.064139, peaks[0][0], places=5)


if __name__ == "__main__":
    unittest.main()
