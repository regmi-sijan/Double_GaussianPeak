"""Linear plus exponential sideband subtraction. Adapted from the archived notebook."""
from .common import load_histogram, options, export_result
import ROOT

def main(argv=None):
    args = options("sideband", argv)
    ROOT.gROOT.SetBatch(True)
    hist = load_histogram(args.input, args.histogram)
    import math
    from array import array

    # ------------------------------------------------------------
    # Settings
    # ------------------------------------------------------------
    signal_min = 820.0
    signal_max = 1100.0

    number_of_left_bins = 50
    number_of_right_bins = 50

    # Centering x improves numerical stability
    x_reference = 0.5 * (signal_min + signal_max)

    # ------------------------------------------------------------
    # Find sideband bins
    # ------------------------------------------------------------
    nbins = hist.GetNbinsX()

    left_candidates = [
        bin_number
        for bin_number in range(1, nbins + 1)
        if hist.GetBinCenter(bin_number) < signal_min
    ]

    right_candidates = [
        bin_number
        for bin_number in range(1, nbins + 1)
        if hist.GetBinCenter(bin_number) > signal_max
    ]

    if len(left_candidates) < number_of_left_bins:
        raise ValueError(
            f"Not enough bins below {signal_min}."
        )

    if len(right_candidates) < number_of_right_bins:
        raise ValueError(
            f"Not enough bins above {signal_max}."
        )

    # Nearest selected bins below the signal interval
    left_bins = left_candidates[-number_of_left_bins:]

    # Nearest selected bins above the signal interval
    right_bins = right_candidates[:number_of_right_bins]

    sideband_bins = left_bins + right_bins

    print("Left sideband:")
    print(
        f"  Number of bins = {len(left_bins)}\n"
        f"  Range = {hist.GetBinCenter(left_bins[0]):.2f} "
        f"to {hist.GetBinCenter(left_bins[-1]):.2f}"
    )

    print("\nRight sideband:")
    print(
        f"  Number of bins = {len(right_bins)}\n"
        f"  Range = {hist.GetBinCenter(right_bins[0]):.2f} "
        f"to {hist.GetBinCenter(right_bins[-1]):.2f}"
    )

    # ------------------------------------------------------------
    # Create a graph containing only the sideband bins
    # ------------------------------------------------------------
    x_values = []
    y_values = []
    x_errors = []
    y_errors = []

    for bin_number in sideband_bins:
        x = hist.GetBinCenter(bin_number)
        y = hist.GetBinContent(bin_number)

        x_error = 0.0  # Bin centers are known; widths are not measurement errors.
        y_error = hist.GetBinError(bin_number)

        if y_error <= 0:
            y_error = math.sqrt(max(abs(y), 1.0))

        x_values.append(x)
        y_values.append(y)
        x_errors.append(x_error)
        y_errors.append(y_error)

    sideband_graph = ROOT.TGraphErrors(
        len(x_values),
        array("d", x_values),
        array("d", y_values),
        array("d", x_errors),
        array("d", y_errors)
    )

    sideband_graph.SetName("background_sideband_graph_pol1_exp")
    sideband_graph.SetTitle("Background sidebands")
    sideband_graph.SetMarkerStyle(20)
    sideband_graph.SetMarkerSize(0.9)
    sideband_graph.SetMarkerColor(ROOT.kBlack)
    sideband_graph.SetLineColor(ROOT.kBlack)

    # ------------------------------------------------------------
    # Background model: pol1 + exponential
    #
    # B(x) = p0 + p1*dx + A*exp(k*dx)
    #
    # dx = x - x_reference
    #
    # [0] = polynomial intercept p0
    # [1] = polynomial slope p1
    # [2] = exponential amplitude A
    # [3] = exponential slope k
    # ------------------------------------------------------------
    def background_model(x, par):
        dx = x[0] - x_reference

        linear_part = par[0] + par[1] * dx
        exponential_part = par[2] * math.exp(par[3] * dx)

        return linear_part + exponential_part


    background_fit_min = hist.GetBinLowEdge(left_bins[0])

    background_fit_max = (
        hist.GetBinLowEdge(right_bins[-1])
        + hist.GetBinWidth(right_bins[-1])
    )

    background_fit = ROOT.TF1(
        "pol1_plus_exponential_background",
        background_model,
        background_fit_min,
        background_fit_max,
        4
    )

    background_fit.SetParName(0, "Polynomial p0")
    background_fit.SetParName(1, "Polynomial p1")
    background_fit.SetParName(2, "Exponential amplitude")
    background_fit.SetParName(3, "Exponential slope")

    # ------------------------------------------------------------
    # Initial estimates
    # ------------------------------------------------------------
    left_y_values = [
        hist.GetBinContent(bin_number)
        for bin_number in left_bins
    ]

    right_y_values = [
        hist.GetBinContent(bin_number)
        for bin_number in right_bins
    ]

    left_average = sum(left_y_values) / len(left_y_values)
    right_average = sum(right_y_values) / len(right_y_values)

    exponential_amplitude_guess = max(
        right_average - left_average,
        1.0
    )

    initial_parameters = [
        max(left_average, 0.0),        # p0
        0.0,                           # p1
        exponential_amplitude_guess,   # exponential amplitude
        0.005                          # exponential slope
    ]

    for index, value in enumerate(initial_parameters):
        background_fit.SetParameter(index, value)

    # Exponential amplitude must be nonnegative
    background_fit.SetParLimits(
        2,
        0.0,
        max(10.0 * max(y_values), 1000.0)
    )

    # Limit exponential slope to prevent numerical divergence
    background_fit.SetParLimits(
        3,
        -0.05,
        0.05
    )

    # ------------------------------------------------------------
    # Fit the 100 selected sideband bins
    # ------------------------------------------------------------
    background_fit_result = sideband_graph.Fit(
        background_fit,
        "RS"
    )

    # ------------------------------------------------------------
    # Select bins in the signal interval
    # ------------------------------------------------------------
    signal_bins = [
        bin_number
        for bin_number in range(1, nbins + 1)
        if (
            hist.GetBinCenter(bin_number) >= signal_min
            and hist.GetBinCenter(bin_number) <= signal_max
        )
    ]

    if not signal_bins:
        raise ValueError(
            f"No histogram bins found between "
            f"{signal_min} and {signal_max}."
        )

    first_signal_bin = signal_bins[0]
    last_signal_bin = signal_bins[-1]

    # ------------------------------------------------------------
    # Create background and subtracted-signal histograms
    # ------------------------------------------------------------
    background_hist = hist.Clone(
        "estimated_pol1_exp_background_hist"
    )
    background_hist.Reset("ICES")

    signal_hist = hist.Clone(
        "pol1_exp_background_subtracted_signal"
    )
    signal_hist.Reset("ICES")

    observed_yield = 0.0
    background_yield = 0.0
    signal_yield = 0.0
    data_variance = 0.0

    for bin_number in signal_bins:
        x = hist.GetBinCenter(bin_number)

        observed_content = hist.GetBinContent(bin_number)
        observed_error = hist.GetBinError(bin_number)

        estimated_background = background_fit.Eval(x)

        background_subtracted_content = (
            observed_content - estimated_background
        )

        background_hist.SetBinContent(
            bin_number,
            estimated_background
        )

        signal_hist.SetBinContent(
            bin_number,
            background_subtracted_content
        )

        # Histogram uncertainty only
        signal_hist.SetBinError(
            bin_number,
            observed_error
        )

        observed_yield += observed_content
        background_yield += estimated_background
        signal_yield += background_subtracted_content
        data_variance += observed_error**2

    data_only_signal_error = math.sqrt(data_variance)

    # ------------------------------------------------------------
    # Fit quality
    # ------------------------------------------------------------
    chi_square = background_fit.GetChisquare()
    ndf = background_fit.GetNDF()

    reduced_chi_square = (
        chi_square / ndf
        if ndf > 0
        else float("nan")
    )

    # ------------------------------------------------------------
    # Draw original histogram and background fit
    # ------------------------------------------------------------
    canvas_background = ROOT.TCanvas(
        "canvas_pol1_exp_background",
        "pol1 plus exponential background",
        600,
        400
    )

    hist.SetLineColor(ROOT.kBlue + 1)
    hist.SetLineWidth(2)
    hist.SetFillColorAlpha(ROOT.kAzure - 9, 0.25)

    display_min = (
        hist.GetBinLowEdge(left_bins[0]) - 20
    )

    display_max = (
        hist.GetBinLowEdge(right_bins[-1])
        + hist.GetBinWidth(right_bins[-1])
        + 20
    )

    hist.GetXaxis().SetRangeUser(
        display_min,
        display_max
    )

    hist.GetXaxis().SetTitle("Channel")
    hist.GetYaxis().SetTitle("Counts")
    hist.Draw("HIST")

    background_fit.SetLineColor(ROOT.kRed)
    background_fit.SetLineWidth(3)
    background_fit.Draw("SAME")

    sideband_graph.Draw("P SAME")

    # Boundaries of excluded region
    left_boundary = ROOT.TLine(
        signal_min,
        0,
        signal_min,
        hist.GetMaximum()
    )

    right_boundary = ROOT.TLine(
        signal_max,
        0,
        signal_max,
        hist.GetMaximum()
    )

    for line in [left_boundary, right_boundary]:
        line.SetLineColor(ROOT.kMagenta + 1)
        line.SetLineStyle(7)
        line.SetLineWidth(2)
        line.Draw("SAME")

    legend_background = ROOT.TLegend(
        0.50,
        0.67,
        0.88,
        0.87
    )

    legend_background.SetBorderSize(0)
    legend_background.SetFillStyle(0)

    legend_background.AddEntry(
        hist,
        "Original histogram",
        "lf"
    )

    legend_background.AddEntry(
        sideband_graph,
        "50 left + 50 right bins",
        "p"
    )

    legend_background.AddEntry(
        background_fit,
        "pol1 + exponential",
        "l"
    )

    legend_background.AddEntry(
        left_boundary,
        "Excluded region: 820-1100",
        "l"
    )

    legend_background.Draw()
    canvas_background.Draw()

    # ------------------------------------------------------------
    # Draw the background-subtracted signal
    # ------------------------------------------------------------
    canvas_signal = ROOT.TCanvas(
        "canvas_pol1_exp_subtracted_signal",
        "Background-subtracted signal",
        600,
        400
    )

    signal_hist.SetLineColor(ROOT.kBlue + 1)
    signal_hist.SetLineWidth(2)
    signal_hist.SetFillColorAlpha(ROOT.kAzure - 9, 0.35)

    signal_hist.GetXaxis().SetRangeUser(
        signal_min,
        signal_max
    )

    signal_hist.GetXaxis().SetTitle("Channel")
    signal_hist.GetYaxis().SetTitle(
        "Background-subtracted counts"
    )

    signal_hist.Draw("HIST E")

    zero_line = ROOT.TLine(
        signal_min,
        0,
        signal_max,
        0
    )

    zero_line.SetLineColor(ROOT.kBlack)
    zero_line.SetLineStyle(2)
    zero_line.Draw("SAME")

    canvas_signal.Draw()

    # ------------------------------------------------------------
    # Print results
    # ------------------------------------------------------------
    print("\n========== BACKGROUND FIT RESULTS ==========")

    print("\nBackground parameters:")
    print(
        f"  Polynomial p0 = "
        f"{background_fit.GetParameter(0):.6g} "
        f"± {background_fit.GetParError(0):.6g}"
    )

    print(
        f"  Polynomial p1 = "
        f"{background_fit.GetParameter(1):.6g} "
        f"± {background_fit.GetParError(1):.6g}"
    )

    print(
        f"  Exponential amplitude = "
        f"{background_fit.GetParameter(2):.6g} "
        f"± {background_fit.GetParError(2):.6g}"
    )

    print(
        f"  Exponential slope = "
        f"{background_fit.GetParameter(3):.6g} "
        f"± {background_fit.GetParError(3):.6g}"
    )

    print("\nBackground fit quality:")
    print(f"  Chi-square     = {chi_square:.3f}")
    print(f"  NDF            = {ndf}")
    print(f"  Chi-square/NDF = {reduced_chi_square:.3f}")
    print(f"  Fit status     = {int(background_fit_result)}")

    print(f"\n========== YIELD FROM {signal_min:g} TO {signal_max:g} ==========")

    print(
        f"  First included bin center = "
        f"{hist.GetBinCenter(first_signal_bin):.3f}"
    )

    print(
        f"  Last included bin center  = "
        f"{hist.GetBinCenter(last_signal_bin):.3f}"
    )

    print(
        f"  Number of signal bins     = "
        f"{len(signal_bins)}"
    )

    print(
        f"\n  Observed histogram yield    = "
        f"{observed_yield:.2f}"
    )

    print(
        f"  Estimated background        = "
        f"{background_yield:.2f}"
    )

    print(
        f"  Background-subtracted yield = "
        f"{signal_yield:.2f}"
    )

    print(
        f"  Data-only statistical error = "
        f"{data_only_signal_error:.2f}"
    )

    if int(background_fit_result) != 0:
        print(
            "\nWarning: the background fit did not converge."
        )

    export_result(args, hist, background_fit, background_fit_result, [canvas_background, canvas_signal], {
        "observed_yield": observed_yield,
        "background_yield": background_yield,
        "signal_yield": signal_yield,
        "data_only_signal_error": data_only_signal_error,
        "signal_interval_bin_centers": [signal_min, signal_max],
        "uncertainty_note": "Background covariance and model uncertainty are not included.",
    })

if __name__ == "__main__":
    main()
