"""Regression and failure-path checks for the supplied spectrum."""
import json
import subprocess
import sys
import tempfile
import unittest
from array import array
from pathlib import Path
import ROOT
from analysis.common import REPO, load_histogram


class AnalysisTests(unittest.TestCase):
    def test_supplied_histogram(self):
        hist = load_histogram(REPO / 'data/mvmelst_008_e10.root', 'h_e10')
        self.assertEqual(hist.GetNbinsX(), 8192)
        self.assertEqual(hist.Integral(), 943707)
        self.assertAlmostEqual(hist.GetBinWidth(1), 8191/8192)

    def test_missing_input_and_object(self):
        with self.assertRaises(FileNotFoundError):
            load_histogram(REPO / 'does-not-exist.root', 'h_e10')
        with self.assertRaises(ValueError):
            load_histogram(REPO / 'data/mvmelst_008_e10.root', 'missing')

    def test_variable_bins_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'variable.root'
            f = ROOT.TFile(str(path), 'RECREATE')
            hist = ROOT.TH1D('variable', '', 2, array('d', [0, 1, 3]))
            hist.Write()
            f.Close()
            with self.assertRaises(ValueError):
                load_histogram(path, 'variable')

    def test_fits_and_strict_diagnostics(self):
        with tempfile.TemporaryDirectory() as tmp:
            for method in ['three_gaussian', 'crystal_ball', 'sideband']:
                with self.subTest(method=method):
                    output = Path(tmp) / method
                    run = subprocess.run([sys.executable, '-m', f'analysis.{method}',
                                          '--output', str(output), '--strict'],
                                         cwd=REPO, capture_output=True, text=True)
                    self.assertTrue((output / 'fit.json').exists(), run.stdout + run.stderr)
                    result = json.loads((output / 'fit.json').read_text())
                    self.assertTrue((output / 'figure_1.png').exists())
                    if method == 'three_gaussian':
                        self.assertEqual(run.returncode, 0, run.stderr)
                        self.assertTrue(result['passes_numerical_checks'])
                        self.assertAlmostEqual(result['yields']['small_peak_count_3sigma'], 21902.69, delta=5)
                    else:
                        self.assertNotEqual(run.returncode, 0)
                        self.assertFalse(result['passes_numerical_checks'])
                        self.assertTrue(result['warnings'])


if __name__ == '__main__':
    unittest.main()
