import tempfile
import unittest
from pathlib import Path

import numpy as np
from matchms import Spectrum
from matchms.similarity import CosineGreedy, ModifiedCosine

from common import read_json, write_json
from preprocess import preprocess


class AnalysisTests(unittest.TestCase):
    def test_shifted_peaks(self):
        a = Spectrum(mz=np.array([50., 70., 90.]), intensities=np.array([1., .5, .2]),
                     metadata={"precursor_mz": 120.}, metadata_harmonization=False)
        b = Spectrum(mz=np.array([60., 80., 100.]), intensities=np.array([1., .5, .2]),
                     metadata={"precursor_mz": 130.}, metadata_harmonization=False)
        self.assertEqual(float(CosineGreedy(tolerance=.02).pair(a, b)["score"]), 0.)
        result = ModifiedCosine(tolerance=.02).pair(a, b)
        self.assertAlmostEqual(float(result["score"]), 1.)
        self.assertEqual(int(result["matches"]), 3)
        self.assertAlmostEqual(float(CosineGreedy(tolerance=.02).pair(a, a)["score"]), 1.)

    def test_tampered_input(self):
        root = Path(__file__).resolve().parents[1]
        manifest = read_json(root / "data/manifest.json")
        manifest["records"][0]["sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            write_json(path / "manifest.json", manifest)
            with self.assertRaisesRegex(ValueError, "Checksum mismatch"):
                preprocess(root / "data/raw", path / "manifest.json", read_json(root / "params.json"), path / "out")

    def test_pair_integrity(self):
        import pandas as pd
        root = Path(__file__).resolve().parents[1]
        records = read_json(root / "expected/preprocess/spectra.json")
        table = pd.read_csv(root / "expected/report/pairs.csv")
        self.assertEqual(len(table), len(records) * (len(records) - 1) // 2)
        self.assertEqual(len({r["canonical_smiles"] for r in records}), len(records))
        self.assertFalse(table.duplicated(["id_a", "id_b"]).any())
        self.assertFalse((table.id_a == table.id_b).any())
        self.assertTrue(np.isfinite(table[["cosine", "modified_cosine", "tanimoto"]]).all().all())
        self.assertTrue(table[["cosine", "modified_cosine", "tanimoto"]].apply(lambda s: s.between(0, 1 + 1e-9).all()).all())
        for r in records:
            peaks = np.array(r["peaks"])
            self.assertTrue((np.diff(peaks[:, 0]) > 0).all())
            self.assertTrue((peaks[:, 0] < r["precursor_mz"] - 1.5).all())
            self.assertGreaterEqual(len(peaks), 5)
            self.assertAlmostEqual(peaks[:, 1].max(), 1.)


if __name__ == "__main__":
    unittest.main()
