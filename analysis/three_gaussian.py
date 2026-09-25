"""Three Gaussians and quadratic background. Adapted from the archived notebook."""
from .common import load_histogram, options, export_result
import ROOT

def main(argv=None):
    args = options("three_gaussian", argv)
    ROOT.gROOT.SetBatch(True)
    hist = load_histogram(args.input, args.histogram)


    # ------------------------------------------------------------
    # Fit settings
    # ------------------------------------------------------------
    fit_min = 750
    fit_max = 1600

    # ------------------------------------------------------------
    # Model:
    #   Small Gaussian + Gaussian 1 + Gaussian 2 + pol2 background
    #
    # [0], [1], [2]    = small peak: amplitude, mean, sigma
    # [3], [4], [5]    = peak 1: amplitude, mean, sigma
    # [6], [7], [8]    = peak 2: amplitude, mean, sigma
    # [9], [10], [11]  = background: p0, p1, p2
    # ------------------------------------------------------------
    fit_func = ROOT.TF1(
        "three_gaussian_pol2_fit",
        "gaus(0) + gaus(3) + gaus(6) + pol2(9)",
        fit_min,
        fit_max
    )

    # Initial parameter estimates
    initial_parameters = [
        200, 980, 55,       # Small peak
        2000, 1240, 45,     # Large peak 1
        2300, 1360, 65,     # Large peak 2
        10, 0, 0            # Quadratic background
    ]

    # Set parameters individually to avoid PyROOT's argument limit
    for index, value in enumerate(initial_parameters):
        fit_func.SetParameter(index, value)

    # Give the parameters descriptive names
    parameter_names = [
        "Small amplitude",
        "Small mean",
        "Small sigma",
        "Peak 1 amplitude",
        "Peak 1 mean",
        "Peak 1 sigma",
        "Peak 2 amplitude",
        "Peak 2 mean",
        "Peak 2 sigma",
        "Background p0",
        "Background p1",
        "Background p2"
    ]

    for index, name in enumerate(parameter_names):
        fit_func.SetParName(index, name)

    # ------------------------------------------------------------
    # Parameter constraints
    # ------------------------------------------------------------

    # Small peak
    fit_func.SetParLimits(0, 0, 1000)
    fit_func.SetParLimits(1, 900, 1050)
    fit_func.SetParLimits(2, 10, 120)

    # Large peak 1
    fit_func.SetParLimits(3, 0, 5000)
    fit_func.SetParLimits(4, 1180, 1300)
    fit_func.SetParLimits(5, 5, 150)

    # Large peak 2
    fit_func.SetParLimits(6, 0, 5000)
    fit_func.SetParLimits(7, 1300, 1450)
    fit_func.SetParLimits(8, 5, 150)

    # ------------------------------------------------------------
    # Perform the fit
    #
    # R = use the specified fit range
    # S = return the fit-result object
    # ------------------------------------------------------------
    fit_result = hist.Fit(fit_func, "RS")

    # ------------------------------------------------------------
    # Create individual fitted components
    # ------------------------------------------------------------

    # Small Gaussian
    small_gaussian = ROOT.TF1(
        "small_gaussian_component",
        "gaus",
        fit_min,
        fit_max
    )

    small_gaussian.SetParameters(
        fit_func.GetParameter(0),
        fit_func.GetParameter(1),
        fit_func.GetParameter(2)
    )

    small_gaussian.SetLineColor(ROOT.kOrange + 7)
    small_gaussian.SetLineStyle(2)
    small_gaussian.SetLineWidth(3)

    # Large Gaussian 1
    gaussian1 = ROOT.TF1(
        "large_gaussian1_component",
        "gaus",
        fit_min,
        fit_max
    )

    gaussian1.SetParameters(
        fit_func.GetParameter(3),
        fit_func.GetParameter(4),
        fit_func.GetParameter(5)
    )

    gaussian1.SetLineColor(ROOT.kGreen + 2)
    gaussian1.SetLineStyle(2)
    gaussian1.SetLineWidth(3)

    # Large Gaussian 2
    gaussian2 = ROOT.TF1(
        "large_gaussian2_component",
        "gaus",
        fit_min,
        fit_max
    )

    gaussian2.SetParameters(
        fit_func.GetParameter(6),
        fit_func.GetParameter(7),
        fit_func.GetParameter(8)
    )

    gaussian2.SetLineColor(ROOT.kMagenta + 1)
    gaussian2.SetLineStyle(2)
    gaussian2.SetLineWidth(3)

    # Quadratic background
    background = ROOT.TF1(
        "pol2_background_component",
        "pol2",
        fit_min,
        fit_max
    )

    background.SetParameters(
        fit_func.GetParameter(9),
        fit_func.GetParameter(10),
        fit_func.GetParameter(11)
    )

    background.SetLineColor(ROOT.kGray + 2)
    background.SetLineStyle(3)
    background.SetLineWidth(3)

    # ------------------------------------------------------------
    # Draw histogram and fitted components
    # ------------------------------------------------------------
    canvas_fit = ROOT.TCanvas(
        "canvas_three_gaussian_pol2",
        "Three Gaussians with quadratic background",
        700,
        500
    )

    hist.SetLineColor(ROOT.kBlue + 1)
    hist.SetLineWidth(2)
    hist.SetFillColorAlpha(ROOT.kAzure - 9, 0.30)

    hist.GetXaxis().SetRangeUser(600, 2000)
    hist.GetXaxis().SetTitle("Channel")
    hist.GetYaxis().SetTitle("Counts")
    hist.Draw("HIST")

    # Total fitted model
    fit_func.SetLineColor(ROOT.kRed)
    fit_func.SetLineWidth(3)
    fit_func.Draw("SAME")

    # Individual fitted components
    small_gaussian.Draw("SAME")
    gaussian1.Draw("SAME")
    gaussian2.Draw("SAME")
    background.Draw("SAME")

    # Legend
    legend = ROOT.TLegend(0.58, 0.57, 0.88, 0.87)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.AddEntry(hist, "Data", "lf")
    legend.AddEntry(fit_func, "Total fit", "l")
    legend.AddEntry(small_gaussian, "Small Gaussian", "l")
    legend.AddEntry(gaussian1, "Gaussian 1", "l")
    legend.AddEntry(gaussian2, "Gaussian 2", "l")
    legend.AddEntry(background, "Quadratic background", "l")

    legend.Draw()
    canvas_fit.Draw()

    # ------------------------------------------------------------
    # Extract fitted parameters
    # ------------------------------------------------------------

    # Small peak
    small_amplitude = fit_func.GetParameter(0)
    small_amplitude_error = fit_func.GetParError(0)

    small_mean = fit_func.GetParameter(1)
    small_mean_error = fit_func.GetParError(1)

    small_sigma = fit_func.GetParameter(2)
    small_sigma_error = fit_func.GetParError(2)


    # Fit quality
    chi_square = fit_func.GetChisquare()
    ndf = fit_func.GetNDF()

    reduced_chi_square = (
        chi_square / ndf
        if ndf > 0
        else float("nan")
    )

    # ------------------------------------------------------------
    # Print results for dedicated peak
    # ------------------------------------------------------------

    print("\nSmall Gaussian:")
    print(
        f"  Mean      = {small_mean:.2f} "
        f"± {small_mean_error:.2f}"
    )
    print(
        f"  Sigma     = {small_sigma:.2f} "
        f"± {small_sigma_error:.2f}"
    )

    print("\nFit quality:")
    print(f"  Chi-square     = {chi_square:.2f}")
    print(f"  NDF            = {ndf}")
    print(f"  Chi-square/NDF = {reduced_chi_square:.2f}")
    print(f"  Fit status     = {int(fit_result)}")

    # Warn if the fit did not converge
    if int(fit_result) != 0:
        print(
            "\nWarning: the fit did not converge successfully. "
            "Try adjusting the initial estimates or parameter limits."
        )

    import math

    # ------------------------------------------------------------
    # Count events in the small peak within mean ± 3 sigma
    # ------------------------------------------------------------
    small_lower_limit = small_mean - 3.0 * small_sigma
    small_upper_limit = small_mean + 3.0 * small_sigma

    # Histogram bin width, assuming uniform binning
    mean_bin = hist.GetXaxis().FindBin(small_mean)
    bin_width = hist.GetXaxis().GetBinWidth(mean_bin)

    # Integral of the fitted Gaussian between mean - 3 sigma
    # and mean + 3 sigma.
    #
    # Division by bin width converts the TF1 area into histogram counts.
    small_peak_count = (
        small_amplitude
        * small_sigma
        * math.sqrt(2.0 * math.pi)
        * math.erf(3.0 / math.sqrt(2.0))
        / bin_width
    )

    # Total fitted Gaussian count over the full range (-infinity, +infinity)
    small_peak_total_count = (
        small_amplitude
        * small_sigma
        * math.sqrt(2.0 * math.pi)
        / bin_width
    )

    # Fraction of a Gaussian inside ±3 sigma
    fraction_inside_3sigma = math.erf(3.0 / math.sqrt(2.0))

    # Observed histogram counts in the same interval.
    # This includes signal and background.
    first_bin = hist.GetXaxis().FindBin(small_lower_limit)
    last_bin = hist.GetXaxis().FindBin(small_upper_limit)

    observed_count = hist.Integral(first_bin, last_bin)


    print("\nSmall-peak counts:")
    print(
        f"  Integration interval = "
        f"[{small_lower_limit:.2f}, {small_upper_limit:.2f}]"
    )
    print(f"  Histogram bin width  = {bin_width:.4g}")
    print(f"  Gaussian count ±3σ   = {small_peak_count:.2f}")
    print(f"  Full Gaussian count  = {small_peak_total_count:.2f}")
    print(f"  Observed bin count   = {observed_count:.2f}")
    print(
        f"  Gaussian fraction inside ±3σ = "
        f"{100.0 * fraction_inside_3sigma:.3f}%"
    )


    export_result(args, hist, fit_func, fit_result, [canvas_fit], {
        "small_peak_count_3sigma": small_peak_count,
        "small_peak_count_full": small_peak_total_count,
        "observed_count_whole_bins": observed_count,
        "integration_interval": [small_lower_limit, small_upper_limit],
    })

if __name__ == "__main__":
    main()
