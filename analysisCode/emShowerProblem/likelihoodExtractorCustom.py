"""
likelihoodExtractorCustom.py
-----------------------------
PyROOT rewrite of likelihoodExtractor.py: parses and plots the photon ID
likelihood PDFs stored in PandoraLikelihoodData12EBin.xml (from
madbaron/SteeringMacros), styled to match profileCompCustom.jl /
rzCrossSectionCustom.jl (transparent background, custom color scheme,
.svg output instead of matplotlib PNGs).

Usage
-----
  python likelihoodExtractorCustom.py [path/to/PandoraLikelihoodData12EBin.xml]
                                      [--output-dir likelihoodPlotsCustom]

Output files (one canvas each, saved as SVG)
------------
  pandora_<variable>_all_ebins.svg   - one canvas per variable, all energy bins
  pandora_all_vars_ebin<N>.svg       - one canvas per energy bin, all variables
  pandora_likelihood_ratios.svg      - log(S/B) ratio for every variable x energy bin
"""

import sys
import os
import argparse
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("xml_path", nargs="?", default="PandoraLikelihoodData12EBin.xml")
parser.add_argument("--output-dir", default="likelihoodPlotsCustom")
args = parser.parse_args()

if not os.path.exists(args.xml_path):
    print(f"ERROR: Cannot find {args.xml_path}")
    print("Usage: python likelihoodExtractorCustom.py [path/to/PandoraLikelihoodData12EBin.xml]")
    sys.exit(1)

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# --- Color scheme (matches profileCompCustom.jl / rzCrossSectionCustom.jl) ---
COLORS_HEX = {
    "purple":  "#101c4d",
    "nominal": "#e38153",
    "vacuum":  "#d14382",
    "blue":    "#0052c5",
    "gold":    "#e3ab0f",
    "yellow":  "#ffe436",
    "text":    "#eeeeee",
    "axis":    "#eeeeee",
}


def hex_to_rgb01(hexstr):
    hexstr = hexstr.lstrip("#")
    return (int(hexstr[0:2], 16) / 255.0,
            int(hexstr[2:4], 16) / 255.0,
            int(hexstr[4:6], 16) / 255.0)


def hex_to_root_color(hexstr):
    r, g, b = hex_to_rgb01(hexstr)
    return ROOT.TColor.GetColor(r, g, b)


COLORS = {k: hex_to_root_color(v) for k, v in COLORS_HEX.items()}
TEXT_COLOR = COLORS["text"]
AXIS_COLOR = COLORS["axis"]

# GRADIENT: purple -> blue -> vacuum(pink), mirroring the Julia
# cgrad([:purple, :blue, :vacuum]) used for the pseudolayer/energy-bin colorbars.
GRADIENT_STOPS = [hex_to_rgb01(COLORS_HEX["purple"]),
                   hex_to_rgb01(COLORS_HEX["blue"]),
                   hex_to_rgb01(COLORS_HEX["vacuum"])]


def gradient_color(frac):
    """frac in [0, 1] -> ROOT color index, interpolated across GRADIENT_STOPS."""
    frac = min(max(frac, 0.0), 1.0)
    n_segs = len(GRADIENT_STOPS) - 1
    seg = min(int(frac * n_segs), n_segs - 1)
    local = frac * n_segs - seg
    r0, g0, b0 = GRADIENT_STOPS[seg]
    r1, g1, b1 = GRADIENT_STOPS[seg + 1]
    r = r0 + (r1 - r0) * local
    g = g0 + (g1 - g0) * local
    b = b0 + (b1 - b0) * local
    return ROOT.TColor.GetColor(r, g, b)


# ── pretty names and axis labels ─────────────────────────────────────────────

VARIABLE_LABELS = {
    "PeakRms":                ("Peak RMS",             "Transverse RMS at shower peak  [cm]"),
    "LongProfileStart":       ("Long. Profile Start",  "Shower start layer"),
    "LongProfileDiscrepancy": ("Long. Profile Discr.", "Longitudinal profile discrepancy"),
    "PeakEnergyFraction":     ("Peak Energy Fraction", "Fraction of energy in shower peak"),
    "RmsRatio":               ("RMS Ratio",            "Peak RMS / overall cluster RMS"),
    "MinDistanceToTrack":     ("Min Track Distance",   "Minimum distance to nearest track  [mm]"),
}

CANONICAL_ORDER = ["PeakRms", "LongProfileStart", "LongProfileDiscrepancy",
                   "PeakEnergyFraction", "RmsRatio", "MinDistanceToTrack"]


# ── parsing ───────────────────────────────────────────────────────────────────

def load_xml(path):
    with open(path) as f:
        raw = f.read()
    # File has no single root element - wrap it
    return ET.fromstring(f"<root>{raw}</root>")


def parse_global(root):
    n_ebins = int(root.find("NEnergyBins").text)
    e_edges = list(map(float, root.find("EnergyBinLowerEdges").text.split()))
    n_sig   = list(map(int,   root.find("NSignalEvents").text.split()))
    n_bkg   = list(map(int,   root.find("NBackgroundEvents").text.split()))
    return n_ebins, e_edges, n_sig, n_bkg


def parse_histogram(root, tag):
    """Return (nbins, lo, hi, contents) for the named tag, or None."""
    el = root.find(tag)
    if el is None:
        return None
    n = int(el.find("NBinsX").text)
    lo = float(el.find("XLow").text)
    hi = float(el.find("XHigh").text)
    raw = list(map(float, el.find("BinContents").text.split()))
    contents = raw[:n]  # XML sometimes stores n+1 values (overflow sentinel)
    return n, lo, hi, contents


def load_all_histograms(root, variables, n_ebins):
    """data[variable][ebin] = {"sig": (n, lo, hi, contents), "bkg": (...)}"""
    data = {}
    for var in variables:
        data[var] = {}
        for i in range(n_ebins):
            sig = parse_histogram(root, f"PhotonSig{var}_{i}")
            bkg = parse_histogram(root, f"PhotonBkg{var}_{i}")
            if sig and bkg:
                data[var][i] = {"sig": sig, "bkg": bkg}
    return data


def ebin_label(i, e_edges, n_ebins):
    lo = e_edges[i]
    hi = e_edges[i + 1] if i + 1 < len(e_edges) else None
    return f"{lo:g}-{hi:g} GeV" if hi is not None else f">{lo:g} GeV"


def make_hist(name, hist_tuple, color):
    n, lo, hi, contents = hist_tuple
    h = ROOT.TH1D(name, "", n, lo, hi)
    for i, val in enumerate(contents):
        h.SetBinContent(i + 1, val)
    h.SetLineColor(color)
    h.SetLineWidth(2)
    h.SetFillColorAlpha(color, 0.25)
    return h


def style_pad(pad):
    pad.SetFillColorAlpha(ROOT.kWhite, 0)
    pad.SetFillStyle(4000)
    pad.SetFrameFillStyle(4000)
    pad.SetTickx(1)
    pad.SetTicky(1)
    pad.SetLeftMargin(0.16)
    pad.SetBottomMargin(0.16)
    pad.SetTopMargin(0.12)
    pad.SetRightMargin(0.05)


def style_hist_axes(h, xlabel, title_size=0.075, label_size=0.06):
    h.SetTitle("")
    xax, yax = h.GetXaxis(), h.GetYaxis()
    xax.SetTitle(xlabel)
    yax.SetTitle("Density")
    for axis in (xax, yax):
        axis.SetTitleColor(TEXT_COLOR)
        axis.SetLabelColor(TEXT_COLOR)
        axis.SetAxisColor(AXIS_COLOR)
        axis.SetTitleSize(title_size)
        axis.SetLabelSize(label_size)
        axis.SetNdivisions(510)
    xax.SetTitleOffset(1.1)
    yax.SetTitleOffset(1.3)


def draw_suptitle(canvas, text):
    lbl = ROOT.TLatex()
    lbl.SetNDC(True)
    lbl.SetTextColor(TEXT_COLOR)
    lbl.SetTextAlign(21)
    lbl.SetTextSize(0.035)
    lbl.DrawLatex(0.5, 0.985, f"#font[62]{{{text}}}")


# ── plot 1: one canvas per variable, all energy bins ─────────────────────────

def plot_per_variable(data, variables, n_ebins, e_edges, out_dir):
    ncols = 4
    nrows = (n_ebins + ncols - 1) // ncols
    keep = []

    for var in variables:
        nice, xlabel = VARIABLE_LABELS.get(var, (var, var))
        c = ROOT.TCanvas(f"c_{var}", nice, ncols * 380, nrows * 320)
        c.SetFillColorAlpha(ROOT.kWhite, 0)
        c.SetFillStyle(4000)
        c.Divide(ncols, nrows, 0.002, 0.002)

        for i in range(n_ebins):
            pad = c.cd(i + 1)
            style_pad(pad)
            hists = data[var].get(i)
            if not hists:
                continue
            color = gradient_color(i / max(n_ebins - 1, 1))
            h_sig = make_hist(f"h_sig_{var}_{i}", hists["sig"], COLORS["blue"])
            h_bkg = make_hist(f"h_bkg_{var}_{i}", hists["bkg"], COLORS["vacuum"])
            style_hist_axes(h_sig, xlabel)
            ymax = max(h_sig.GetMaximum(), h_bkg.GetMaximum()) * 1.2
            h_sig.SetMaximum(ymax)
            h_sig.SetMinimum(0)
            h_sig.Draw("HIST")
            h_bkg.Draw("HIST SAME")

            title = ROOT.TLatex()
            title.SetNDC(True)
            title.SetTextColor(color)
            title.SetTextAlign(21)
            title.SetTextSize(0.09)
            title.DrawLatex(0.5, 0.93, ebin_label(i, e_edges, n_ebins))

            if i == 0:
                leg = ROOT.TLegend(0.45, 0.72, 0.92, 0.90)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.SetTextColor(TEXT_COLOR)
                leg.SetTextSize(0.075)
                leg.AddEntry(h_sig, "Signal (#gamma)", "f")
                leg.AddEntry(h_bkg, "Background", "f")
                leg.Draw()
                keep.append(leg)
            keep += [h_sig, h_bkg]

        c.cd(0)
        draw_suptitle(c, f"Photon ID PDFs -- {nice}")

        fname = os.path.join(out_dir, f"pandora_{var}_all_ebins.svg")
        c.SaveAs(fname)
        print(f"  Saved {fname}")
        keep.append(c)


# ── plot 2: one canvas per energy bin, all variables ─────────────────────────

def plot_per_ebin(data, variables, n_ebins, e_edges, out_dir):
    nv = len(variables)
    ncols = 3
    nrows = (nv + ncols - 1) // ncols
    keep = []

    for i in range(n_ebins):
        c = ROOT.TCanvas(f"c_ebin{i}", f"ebin {i}", ncols * 420, nrows * 360)
        c.SetFillColorAlpha(ROOT.kWhite, 0)
        c.SetFillStyle(4000)
        c.Divide(ncols, nrows, 0.005, 0.005)

        for j, var in enumerate(variables):
            pad = c.cd(j + 1)
            style_pad(pad)
            nice, xlabel = VARIABLE_LABELS.get(var, (var, var))
            hists = data[var].get(i)
            if not hists:
                continue
            h_sig = make_hist(f"h_sig_e{i}_{var}", hists["sig"], COLORS["blue"])
            h_bkg = make_hist(f"h_bkg_e{i}_{var}", hists["bkg"], COLORS["vacuum"])
            style_hist_axes(h_sig, xlabel, title_size=0.06, label_size=0.05)
            ymax = max(h_sig.GetMaximum(), h_bkg.GetMaximum()) * 1.2
            h_sig.SetMaximum(ymax)
            h_sig.SetMinimum(0)
            h_sig.Draw("HIST")
            h_bkg.Draw("HIST SAME")

            title = ROOT.TLatex()
            title.SetNDC(True)
            title.SetTextColor(TEXT_COLOR)
            title.SetTextAlign(21)
            title.SetTextSize(0.07)
            title.DrawLatex(0.5, 0.93, nice)

            leg = ROOT.TLegend(0.45, 0.72, 0.92, 0.90)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextColor(TEXT_COLOR)
            leg.SetTextSize(0.06)
            leg.AddEntry(h_sig, "Signal (#gamma)", "f")
            leg.AddEntry(h_bkg, "Background", "f")
            leg.Draw()
            keep += [h_sig, h_bkg, leg]

        c.cd(0)
        draw_suptitle(c, f"Photon ID PDFs -- {ebin_label(i, e_edges, n_ebins)}")

        fname = os.path.join(out_dir, f"pandora_all_vars_ebin{i:02d}.svg")
        c.SaveAs(fname)
        print(f"  Saved {fname}")
        keep.append(c)


# ── plot 3: log(S/B) likelihood ratio overview ───────────────────────────────

def plot_likelihood_ratios(data, variables, n_ebins, e_edges, out_dir):
    import array as pyarray

    nv = len(variables)
    ncols = 3
    nrows = (nv + ncols - 1) // ncols

    c = ROOT.TCanvas("c_ratios", "likelihood ratios", ncols * 430, nrows * 370)
    c.SetFillColorAlpha(ROOT.kWhite, 0)
    c.SetFillStyle(4000)
    c.Divide(ncols, nrows, 0.01, 0.01)
    keep = []

    for j, var in enumerate(variables):
        pad = c.cd(j + 1)
        style_pad(pad)
        nice, xlabel = VARIABLE_LABELS.get(var, (var, var))

        mg = ROOT.TMultiGraph()
        for i in range(n_ebins):
            hists = data[var].get(i)
            if not hists:
                continue
            n_sig, lo_sig, hi_sig, sig = hists["sig"]
            n_bkg, lo_bkg, hi_bkg, bkg = hists["bkg"]
            n = min(n_sig, n_bkg)
            width = (hi_sig - lo_sig) / n
            centres = [lo_sig + width * (k + 0.5) for k in range(n)]
            eps = 1e-9
            ratio = [ROOT.TMath.Log((sig[k] + eps) / (bkg[k] + eps)) for k in range(n)]

            color = gradient_color(i / max(n_ebins - 1, 1))
            g = ROOT.TGraph(n, pyarray.array('d', centres), pyarray.array('d', ratio))
            g.SetLineColor(color)
            g.SetLineWidth(2)
            mg.Add(g, "L")
            keep.append(g)

        mg.Draw("A")
        mg.GetXaxis().SetTitle(xlabel)
        mg.GetYaxis().SetTitle("log(S / B)")
        for axis in (mg.GetXaxis(), mg.GetYaxis()):
            axis.SetTitleColor(TEXT_COLOR)
            axis.SetLabelColor(TEXT_COLOR)
            axis.SetAxisColor(AXIS_COLOR)
            axis.SetTitleSize(0.06)
            axis.SetLabelSize(0.05)
        mg.GetXaxis().SetTitleOffset(1.1)
        mg.GetYaxis().SetTitleOffset(1.3)

        zero_line = ROOT.TLine(mg.GetXaxis().GetXmin(), 0, mg.GetXaxis().GetXmax(), 0)
        zero_line.SetLineColor(AXIS_COLOR)
        zero_line.SetLineStyle(2)
        zero_line.Draw()
        keep += [mg, zero_line]

        title = ROOT.TLatex()
        title.SetNDC(True)
        title.SetTextColor(TEXT_COLOR)
        title.SetTextAlign(21)
        title.SetTextSize(0.07)
        title.DrawLatex(0.5, 0.94, nice)
        keep.append(title)

    # Colour-coded energy-bin legend on the last (possibly empty) pad
    leg_pad_idx = nrows * ncols
    pad = c.cd(leg_pad_idx)
    style_pad(pad)
    leg = ROOT.TLegend(0.02, 0.02, 0.98, 0.98)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextColor(TEXT_COLOR)
    leg.SetTextSize(0.06)
    leg.SetHeader("#font[62]{Energy bin}")
    dummy_graphs = []
    for i in range(n_ebins):
        color = gradient_color(i / max(n_ebins - 1, 1))
        g_dummy = ROOT.TGraph()
        g_dummy.SetLineColor(color)
        g_dummy.SetLineWidth(3)
        dummy_graphs.append(g_dummy)
        leg.AddEntry(g_dummy, ebin_label(i, e_edges, n_ebins), "l")
    leg.Draw()
    keep += dummy_graphs + [leg]

    c.cd(0)
    draw_suptitle(c, "Log-likelihood ratio  log(S/B)  per variable & energy bin")

    fname = os.path.join(out_dir, "pandora_likelihood_ratios.svg")
    c.SaveAs(fname)
    print(f"  Saved {fname}")
    keep.append(c)
    return keep  # keep references alive until process exit


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    print(f"Reading {args.xml_path} ...")
    root = load_xml(args.xml_path)
    n_ebins, e_edges, n_sig, n_bkg = parse_global(root)

    print(f"  {n_ebins} energy bins: {e_edges}")
    print(f"  Signal events per bin:     {n_sig}")
    print(f"  Background events per bin: {n_bkg}")

    variables_found = set()
    for child in root:
        tag = child.tag
        for prefix in ("PhotonSig", "PhotonBkg"):
            if tag.startswith(prefix):
                core = tag[len(prefix):]
                varname = "_".join(core.split("_")[:-1])
                if varname:
                    variables_found.add(varname)

    variables = [v for v in CANONICAL_ORDER if v in variables_found]
    variables += sorted(variables_found - set(CANONICAL_ORDER))
    print(f"  Variables found: {variables}")

    data = load_all_histograms(root, variables, n_ebins)

    os.makedirs(args.output_dir, exist_ok=True)
    print("\nGenerating SVGs ...")
    plot_per_variable(data, variables, n_ebins, e_edges, args.output_dir)
    plot_per_ebin(data, variables, n_ebins, e_edges, args.output_dir)
    _keep = plot_likelihood_ratios(data, variables, n_ebins, e_edges, args.output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
