# Scientific review

## Input and statistical interpretation

The supplied file contains a `TH1F` named `h_e10`, titled `mvmelst_008 - Channel 10`, with 8192 uniform bins over [0, 8191]. Bin width is 8191/8192, not exactly one. The stored x-axis title is Channel. The regular-bin content sum is 943707, whereas `GetEntries()` is 8192. Entries must not be described as the number of detected events; the filling history is unknown.

The histogram has no stored Sumw2 array. ROOT's default bin errors are used; whether Poisson counting assumptions accurately describe the acquisition must be confirmed from provenance. Fits retain the notebook's ROOT chi-square objective and center-sampled functions. ROOT excludes zero-error bins in its default histogram chi-square fit. Consider a bin-integrated Poisson likelihood only after validating nonnegative count data and a nonnegative model.

## Baseline: three Gaussians plus a quadratic

The full model is the sum of three amplitude-normalized Gaussian functions plus a quadratic background. The smaller component has bounds separating it from the two large neighboring components. All 12 parameters are fitted over channels 750–1600.

For a Gaussian of height A and width sigma, the signal yield over mean ±3 sigma is

`A * sigma * sqrt(2*pi) * erf(3/sqrt(2)) / bin_width`.

The conversion assumes a uniform-bin count histogram. The full Gaussian area extends beyond the fit interval and is an extrapolation. Observed histogram sums include complete endpoint bins and background; their intervals do not coincide exactly with continuous function integrals.

The local run reproduces the original baseline: mean about 968.169, sigma about 55.246, chi-square/NDF about 1.242, and status 0. This makes it a useful computational baseline, not proof that the decomposition is unique. Yield uncertainties are not implemented; a complete calculation must include amplitude–width covariance and model/range sensitivity. The saved covariance supports future propagation.

## Crystal Ball: preserve as an exploratory model

The fit uses a right-tailed Crystal Ball signal and a continuous background with a linear left branch, quadratic right branch, and free join. It reproduces a failed fit with zero parameter errors and a join at its upper limit. Do not use these errors or yields as validated results.

In the observed solution, mean + alpha*sigma lies beyond the fitted upper endpoint. The data therefore do not sample the power-law branch, which explains why tail parameters can be unidentifiable. A lower reduced chi-square than the global fit is not evidence of superiority: the fit regions and flexibility differ.

The runnable version uses the clipped mean ±3-sigma interval within the fitted domain. The original notebook's additional reporting cell used an unclipped interval; that contradictory duplicate is retained only in the archive. A Crystal Ball core's ±3-sigma interval does not generally contain the Gaussian 99.73% fraction.

## Sidebands: linear plus exponential, not quadratic

The original section heading said pol2, but its implementation is pol1 plus exponential. The runnable method correctly identifies the linear polynomial. It selects the nearest 50 bin centers on each side of channels 820–1100, fits only those 100 points, and subtracts the estimated background in the excluded interval.

The revised graph uses zero x errors: the bin centers are known and half a bin width is not an uncertainty in their coordinates. Finite-bin averaging would require a separate modeling change. This correction changes the yield from the historical 20809.95 to approximately 20778.91 counts. The exponential slope still reaches its upper limit of 0.05, so the result requires range and model sensitivity studies.

The printed 181.47 error is derived only from the signal-region data errors. Background parameters introduce correlated prediction uncertainty across bins. A complete yield variance requires a covariance contraction of the summed background gradient and the fitted covariance, plus assessment of sideband contamination and model choice. The current error is deliberately labeled data-only.

## Changes made for reproducibility

- Preserved both input files byte for byte, including historical notebook output.
- Split shared-state notebook methods into independent Python modules and notebook entry points.
- Removed install commands and Colab-specific paths from runnable analysis.
- Validated one-dimensional, uniform-bin input and detached it before closing the ROOT file.
- Corrected channel labels, stale sideband ranges, and sideband x errors.
- Exported parameters, covariance, status, bounds, checksums, and figures.
- Flagged fits with failed status, unreliable covariance, or parameters within 0.1% of their allowed range boundaries. The boundary rule is a heuristic for review, not a hypothesis test.

## Next scientific steps

1. Supply acquisition details, calibration, and the intended physical quantity.
2. Define a common signal window when comparing methods.
3. Study residuals, fit-window changes, starting-value dependence, and alternate backgrounds.
4. Propagate yield covariance and study bias/coverage using simulated spectra.
5. Simplify or constrain the Crystal Ball model only with a defensible physical basis.

## References

- [ROOT TH1 fitting and histogram documentation](https://root.cern.ch/doc/master/classTH1.html)
- [ROOT TFitResult diagnostics](https://root.cern.ch/doc/master/classTFitResult.html)
