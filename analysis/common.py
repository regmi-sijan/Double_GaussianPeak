"""Shared input validation and machine-readable fit diagnostics."""
import argparse
import ctypes
import hashlib
import json
import math
from pathlib import Path
import ROOT

REPO = Path(__file__).resolve().parents[1]


def options(method, argv=None):
    parser = argparse.ArgumentParser(description=method)
    parser.add_argument('--input', type=Path, default=REPO / 'data/mvmelst_008_e10.root')
    parser.add_argument('--histogram', default='h_e10')
    parser.add_argument('--output', type=Path, default=REPO / 'results' / method)
    parser.add_argument('--strict', action='store_true', help='Fail on invalid covariance, fit failure, or a parameter near its bound.')
    args = parser.parse_args(argv)
    args.method = method
    return args


def load_histogram(path, name):
    if not Path(path).is_file():
        raise FileNotFoundError(path)
    source = ROOT.TFile.Open(str(path), 'READ')
    if not source or source.IsZombie():
        raise OSError(f'Cannot open ROOT file: {path}')
    try:
        obj = source.Get(name)
        if not obj or not obj.InheritsFrom('TH1') or obj.GetDimension() != 1:
            raise ValueError(f'{name!r} must be a one-dimensional histogram')
        hist = obj.Clone(f'{name}_analysis')
        hist.SetDirectory(0)
    finally:
        source.Close()
    width = hist.GetBinWidth(1)
    if any(not math.isclose(hist.GetBinWidth(i), width, rel_tol=1e-9) for i in range(1, hist.GetNbinsX()+1)):
        raise ValueError('These area-to-count calculations require uniform binning')
    hist.SetStats(False)
    return hist


def export_result(args, hist, function, result, canvases, yields):
    parameters, near_bounds = [], []
    for i in range(function.GetNpar()):
        low, high = ctypes.c_double(), ctypes.c_double()
        function.GetParLimits(i, low, high)
        value = function.GetParameter(i)
        bounded = high.value > low.value
        near = bounded and min(value-low.value, high.value-value) <= 0.001*(high.value-low.value)
        if near:
            near_bounds.append(function.GetParName(i))
        parameters.append(dict(name=function.GetParName(i), value=value,
                               error=function.GetParError(i),
                               bounds=[low.value, high.value] if bounded else None))
    status, covariance_status = int(result), result.CovMatrixStatus()
    reliable = status == 0 and result.IsValid() and covariance_status == 3 and not near_bounds
    warnings = []
    if not reliable:
        warnings.append('Fit requires review: do not interpret reported errors as validated uncertainties.')
    if near_bounds:
        warnings.append('Parameters near bounds: ' + ', '.join(near_bounds))
    report = dict(method=args.method, root_version=ROOT.gROOT.GetVersion(),
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  histogram=args.histogram, fit_range=[function.GetXmin(), function.GetXmax()],
                  fit_status=status, valid=bool(result.IsValid()), covariance_status=covariance_status,
                  passes_numerical_checks=reliable, near_bounds=near_bounds,
                  chi_square=function.GetChisquare(), ndf=function.GetNDF(),
                  parameters=parameters, yields=yields, warnings=warnings,
                  covariance=[[result.CovMatrix(i,j) for j in range(function.GetNpar())]
                              for i in range(function.GetNpar())])
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'fit.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    for index, canvas in enumerate(canvases, 1):
        canvas.SaveAs(str(args.output / f'figure_{index}.png'))
    for warning in warnings:
        print('REVIEW:', warning)
    print(f'Results saved to {args.output}')
    if args.strict and not reliable:
        raise RuntimeError('Fit failed strict numerical checks; diagnostics have been saved.')
