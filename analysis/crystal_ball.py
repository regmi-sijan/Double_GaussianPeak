"""Exploratory Crystal Ball and piecewise background. Adapted from the archived notebook."""
from .common import load_histogram, options, export_result
import ROOT

def main(argv=None):
    args = options("crystal_ball", argv)
    ROOT.gROOT.SetBatch(True)
    hist = load_histogram(args.input, args.histogram)
    import math

    # ------------------------------------------------------------
    # Fit settings
    # ------------------------------------------------------------
    fit_min = 750.0
    fit_max = 1120.0

    # ------------------------------------------------------------
    # Right-tailed Crystal Ball signal
    #
    # par[0] = amplitude
    # par[1] = mean
    # par[2] = sigma
    # par[3] = alpha
    # par[4] = tail exponent n
    # ------------------------------------------------------------
    def crystal_ball_right(x, par):
        amplitude = par[0]
        mean = par[1]
        sigma = abs(par[2])
        alpha = abs(par[3])
        tail_n = par[4]

        if sigma <= 0 or alpha <= 0 or tail_n <= 1:
            return 0.0

        t = (x[0] - mean) / sigma

        # Gaussian core
        if t <= alpha:
            return amplitude * math.exp(-0.5 * t * t)

        # Right-side power-law tail
        A = (
            (tail_n / alpha) ** tail_n
            * math.exp(-0.5 * alpha * alpha)
        )

        B = tail_n / alpha - alpha

        return amplitude * A * (B + t) ** (-tail_n)


    # ------------------------------------------------------------
    # Total model:
    #
    # Crystal Ball signal
    # + pol1 background on the left
    # + pol2 background on the right
    #
    # par[5] = background value at join
    # par[6] = left slope
    # par[7] = right slope
    # par[8] = right quadratic coefficient
    # par[9] = free join point
    #
    # The background is continuous at the join point.
    # ------------------------------------------------------------
    def total_model(x, par):
        signal = crystal_ball_right(x, par)

        background_at_join = par[5]
        left_slope = par[6]
        right_slope = par[7]
        right_quadratic = par[8]
        join_point = par[9]

        dx = x[0] - join_point

        if x[0] <= join_point:
            # Left-side pol1
            background = (
                background_at_join
                + left_slope * dx
            )
        else:
            # Right-side pol2
            background = (
                background_at_join
                + right_slope * dx
                + right_quadratic * dx * dx
            )

        return signal + background


    # ------------------------------------------------------------
    # Create total fit function
    # ------------------------------------------------------------
    fit_func = ROOT.TF1(
        "small_peak_crystal_ball_free_join",
        total_model,
        fit_min,
        fit_max,
        10
    )

    # Parameter names
    parameter_names = [
        "Signal amplitude",      # 0
        "Signal mean",           # 1
        "Signal sigma",          # 2
        "Tail alpha",            # 3
        "Tail n",                # 4
        "Background at join",    # 5
        "Left slope",            # 6
        "Right slope",           # 7
        "Right quadratic",       # 8
        "Join point"             # 9
    ]

    for index, name in enumerate(parameter_names):
        fit_func.SetParName(index, name)

    # ------------------------------------------------------------
    # Initial parameter estimates
    # ------------------------------------------------------------
    initial_parameters = [
        200.0,    # Signal amplitude
        980.0,    # Signal mean
        55.0,     # Signal sigma
        1.5,      # Alpha
        3.0,      # Tail exponent n
        10.0,     # Background value at join
        0.0,      # Left slope
        0.0,      # Right slope
        0.0,      # Right quadratic coefficient
        980.0     # Initial join point
    ]

    for index, value in enumerate(initial_parameters):
        fit_func.SetParameter(index, value)

    # ------------------------------------------------------------
    # Parameter constraints
    # ------------------------------------------------------------
    fit_func.SetParLimits(0, 0.0, 1000.0)   # amplitude
    fit_func.SetParLimits(1, 900.0, 1050.0) # mean
    fit_func.SetParLimits(2, 5.0, 120.0)    # sigma
    fit_func.SetParLimits(3, 0.3, 5.0)      # alpha
    fit_func.SetParLimits(4, 1.01, 30.0)    # tail n

    # Require nonnegative background value at the join
    fit_func.SetParLimits(5, 0.0, 1000.0)

    # Allow the fit to determine the join point within this range
    fit_func.SetParLimits(9, 900.0, 1060.0)

    # ------------------------------------------------------------
    # Perform fit
    #
    # R = use specified fit range
    # S = return fit-result object
    # ------------------------------------------------------------
    fit_result = hist.Fit(fit_func, "RS")

    # ------------------------------------------------------------
    # Extract fitted parameters
    # ------------------------------------------------------------
    amplitude = fit_func.GetParameter(0)
    amplitude_error = fit_func.GetParError(0)

    mean = fit_func.GetParameter(1)
    mean_error = fit_func.GetParError(1)

    sigma = abs(fit_func.GetParameter(2))
    sigma_error = fit_func.GetParError(2)

    alpha = fit_func.GetParameter(3)
    alpha_error = fit_func.GetParError(3)

    tail_n = fit_func.GetParameter(4)
    tail_n_error = fit_func.GetParError(4)

    background_at_join = fit_func.GetParameter(5)
    background_at_join_error = fit_func.GetParError(5)

    left_slope = fit_func.GetParameter(6)
    left_slope_error = fit_func.GetParError(6)

    right_slope = fit_func.GetParameter(7)
    right_slope_error = fit_func.GetParError(7)

    right_quadratic = fit_func.GetParameter(8)
    right_quadratic_error = fit_func.GetParError(8)

    fitted_join_point = fit_func.GetParameter(9)
    fitted_join_error = fit_func.GetParError(9)

    # ------------------------------------------------------------
    # Create the signal-only component
    # ------------------------------------------------------------
    signal_func = ROOT.TF1(
        "small_peak_crystal_ball_signal_free_join",
        crystal_ball_right,
        fit_min,
        fit_max,
        5
    )

    for index in range(5):
        signal_func.SetParameter(
            index,
            fit_func.GetParameter(index)
        )

    signal_func.SetLineColor(ROOT.kMagenta + 1)
    signal_func.SetLineStyle(2)
    signal_func.SetLineWidth(3)

    # ------------------------------------------------------------
    # Create the fitted background-only component
    # ------------------------------------------------------------
    def piecewise_background(x, par):
        background_at_join = par[0]
        left_slope = par[1]
        right_slope = par[2]
        right_quadratic = par[3]
        join_point = par[4]

        dx = x[0] - join_point

        if x[0] <= join_point:
            return background_at_join + left_slope * dx

        return (
            background_at_join
            + right_slope * dx
            + right_quadratic * dx * dx
        )


    background_func = ROOT.TF1(
        "small_peak_background_free_join",
        piecewise_background,
        fit_min,
        fit_max,
        5
    )

    background_func.SetParameters(
        background_at_join,
        left_slope,
        right_slope,
        right_quadratic,
        fitted_join_point
    )

    background_func.SetLineColor(ROOT.kGreen + 2)
    background_func.SetLineStyle(3)
    background_func.SetLineWidth(3)

    # ------------------------------------------------------------
    # Signal integral within fitted mean ± 3 sigma
    # ------------------------------------------------------------
    lower_limit = mean - 3.0 * sigma
    upper_limit = mean + 3.0 * sigma

    # Keep integration within the fitted interval
    integration_min = max(lower_limit, fit_min)
    integration_max = min(upper_limit, fit_max)

    # This conversion assumes uniform histogram bin widths
    mean_bin = hist.GetXaxis().FindBin(mean)
    bin_width = hist.GetXaxis().GetBinWidth(mean_bin)

    # Crystal Ball signal only: polynomial background is excluded
    signal_area_3sigma = signal_func.Integral(
        integration_min,
        integration_max
    )

    signal_count_3sigma = signal_area_3sigma / bin_width

    # Fitted background within the same interval
    background_area_3sigma = background_func.Integral(
        integration_min,
        integration_max
    )

    background_count_3sigma = (
        background_area_3sigma / bin_width
    )

    # Total fitted signal + background
    total_fit_area_3sigma = fit_func.Integral(
        integration_min,
        integration_max
    )

    total_fit_count_3sigma = (
        total_fit_area_3sigma / bin_width
    )

    # Observed histogram count in the same interval
    first_bin = hist.GetXaxis().FindBin(integration_min)
    last_bin = hist.GetXaxis().FindBin(integration_max)

    observed_count_3sigma = hist.Integral(
        first_bin,
        last_bin
    )

    # ------------------------------------------------------------
    # Fit quality
    # ------------------------------------------------------------
    chi_square = fit_func.GetChisquare()
    ndf = fit_func.GetNDF()

    reduced_chi_square = (
        chi_square / ndf
        if ndf > 0
        else float("nan")
    )

    # ------------------------------------------------------------
    # Draw
    # ------------------------------------------------------------
    canvas_small_peak = ROOT.TCanvas(
        "canvas_small_peak_free_join",
        "Small peak with free background join",
        700,
        500
    )

    hist.SetLineColor(ROOT.kBlue + 1)
    hist.SetLineWidth(2)
    hist.SetFillColorAlpha(ROOT.kAzure - 9, 0.30)

    hist.GetXaxis().SetRangeUser(700, 1200)
    hist.GetXaxis().SetTitle("Channel")
    hist.GetYaxis().SetTitle("Counts")
    hist.Draw("HIST")

    # Total model
    fit_func.SetLineColor(ROOT.kRed)
    fit_func.SetLineWidth(3)
    fit_func.Draw("SAME")

    # Signal and background
    signal_func.Draw("SAME")
    background_func.Draw("SAME")

    # Mean ±3 sigma lines
    line_low = ROOT.TLine(
        integration_min,
        0,
        integration_min,
        hist.GetMaximum()
    )

    line_high = ROOT.TLine(
        integration_max,
        0,
        integration_max,
        hist.GetMaximum()
    )

    for line in [line_low, line_high]:
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(7)
        line.SetLineWidth(2)
        line.Draw("SAME")

    # Fitted join-point line
    join_line = ROOT.TLine(
        fitted_join_point,
        0,
        fitted_join_point,
        hist.GetMaximum()
    )

    join_line.SetLineColor(ROOT.kOrange + 7)
    join_line.SetLineStyle(9)
    join_line.SetLineWidth(2)
    join_line.Draw("SAME")

    # Legend
    legend = ROOT.TLegend(0.55, 0.58, 0.88, 0.87)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.AddEntry(hist, "Data", "lf")
    legend.AddEntry(fit_func, "Total fit", "l")
    legend.AddEntry(signal_func, "Crystal Ball signal", "l")
    legend.AddEntry(
        background_func,
        "Piecewise background",
        "l"
    )
    legend.AddEntry(line_low, "Mean #pm 3#sigma", "l")
    legend.AddEntry(join_line, "Fitted join point", "l")

    legend.Draw()
    canvas_small_peak.Draw()

    # ------------------------------------------------------------
    # Print results
    # ------------------------------------------------------------
    print("\n========== SMALL-PEAK FIT RESULTS ==========")

    print("\nCrystal Ball signal:")
    print(
        f"  Amplitude = {amplitude:.3f} "
        f"± {amplitude_error:.3f}"
    )
    print(
        f"  Mean      = {mean:.3f} "
        f"± {mean_error:.3f}"
    )
    print(
        f"  Sigma     = {sigma:.3f} "
        f"± {sigma_error:.3f}"
    )
    print(
        f"  Alpha     = {alpha:.3f} "
        f"± {alpha_error:.3f}"
    )
    print(
        f"  Tail n    = {tail_n:.3f} "
        f"± {tail_n_error:.3f}"
    )

    print("\nPiecewise background:")
    print(
        f"  Fitted join point = {fitted_join_point:.3f} "
        f"± {fitted_join_error:.3f}"
    )
    print(
        f"  Background at join = {background_at_join:.6g} "
        f"± {background_at_join_error:.6g}"
    )
    print(
        f"  Left slope         = {left_slope:.6g} "
        f"± {left_slope_error:.6g}"
    )
    print(
        f"  Right slope        = {right_slope:.6g} "
        f"± {right_slope_error:.6g}"
    )
    print(
        f"  Right quadratic    = {right_quadratic:.6g} "
        f"± {right_quadratic_error:.6g}"
    )

    print("\nCounts within fitted mean ± 3 sigma:")
    print(
        f"  Interval                  = "
        f"[{integration_min:.3f}, {integration_max:.3f}]"
    )
    print(f"  Histogram bin width       = {bin_width:.6g}")
    print(
        f"  Crystal Ball signal count = "
        f"{signal_count_3sigma:.2f}"
    )
    print(
        f"  Fitted background count   = "
        f"{background_count_3sigma:.2f}"
    )
    print(
        f"  Total fitted count        = "
        f"{total_fit_count_3sigma:.2f}"
    )
    print(
        f"  Observed histogram count  = "
        f"{observed_count_3sigma:.2f}"
    )

    print("\nFit quality:")
    print(f"  Chi-square     = {chi_square:.3f}")
    print(f"  NDF            = {ndf}")
    print(f"  Chi-square/NDF = {reduced_chi_square:.3f}")
    print(f"  Fit status     = {int(fit_result)}")

    if int(fit_result) != 0:
        print(
            "\nWarning: the fit did not converge successfully. "
            "Adjust the initial values, limits, or fit range."
        )

    # Check whether the join point reached its allowed boundary
    join_tolerance = 1.0

    if (
        abs(fitted_join_point - 900.0) < join_tolerance
        or abs(fitted_join_point - 1060.0) < join_tolerance
    ):
        print(
            "\nWarning: the fitted join point is near a parameter "
            "boundary. The data may not determine it reliably."
        )

    export_result(args, hist, fit_func, fit_result, [canvas_small_peak], {
        "signal_count_clipped_3sigma": signal_count_3sigma,
        "integration_interval": [integration_min, integration_max],
        "interpretation": "Exploratory; inspect fit validity before using any yield.",
    })

if __name__ == "__main__":
    main()
