"""
averageProfileCustom.py

PyROOT rewrite of the "average profile" canvas from profilePlotter.py,
styled to match profileCompCustom.jl / rzCrossSectionCustom.jl
(transparent background, custom color scheme, .svg output).

Produces a single SVG: <output_dir>/average_profile.svg
  - mean measured vs mean expected longitudinal shower profile,
    averaged over all clusters passing the selection cuts.

Usage
-----
    python averageProfileCustom.py [--input shower_profiles.root]
                                   [--tree ShowerProfiles]
                                   [--output-dir customSVGs]
                                   [--energy-min 0.5]
                                   [--energy-max 1e9]
"""

import argparse
import os
import sys
import array

parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--input",      default="shower_profiles.root")
parser.add_argument("--tree",       default="ShowerProfiles")
parser.add_argument("--output-dir", default="customSVGs")
parser.add_argument("--energy-min", type=float, default=0.5)
parser.add_argument("--energy-max", type=float, default=1e9)
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

BIN_WIDTH = 0.5   # radiation lengths per bin (must match PandoraSettings)

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


def hex_to_root_color(hexstr):
    hexstr = hexstr.lstrip("#")
    r = int(hexstr[0:2], 16) / 255.0
    g = int(hexstr[2:4], 16) / 255.0
    b = int(hexstr[4:6], 16) / 255.0
    return ROOT.TColor.GetColor(r, g, b)


COLORS = {k: hex_to_root_color(v) for k, v in COLORS_HEX.items()}

TEXT_COLOR = COLORS["text"]
AXIS_COLOR = COLORS["axis"]


def make_transparent_canvas(name, title, w, h):
    c = ROOT.TCanvas(name, title, w, h)
    c.SetFillColorAlpha(ROOT.kWhite, 0)
    c.SetFillStyle(4000)
    c.SetFrameFillStyle(4000)
    c.SetLeftMargin(0.13)
    c.SetBottomMargin(0.13)
    c.SetTopMargin(0.10)
    c.SetRightMargin(0.05)
    c.SetTickx(1)
    c.SetTicky(1)
    return c


def style_frame(frame, xlabel, ylabel):
    frame.GetXaxis().SetTitle(xlabel)
    frame.GetYaxis().SetTitle(ylabel)
    for axis in (frame.GetXaxis(), frame.GetYaxis()):
        axis.SetTitleColor(TEXT_COLOR)
        axis.SetLabelColor(TEXT_COLOR)
        axis.SetAxisColor(AXIS_COLOR)
        axis.SetTickLength(0.02)
        axis.SetNdivisions(510)  # 5 major divisions w/ minor ticks, à la IntervalsBetween(5)
    frame.SetTitle("")


def make_graph(x_vals, y_vals, color, line_style=1, line_width=3):
    n = len(x_vals)
    g = ROOT.TGraph(n, array.array('d', x_vals), array.array('d', y_vals))
    g.SetLineColor(color)
    g.SetLineStyle(line_style)
    g.SetLineWidth(line_width)
    return g


# --- Open tree ---
tfile = ROOT.TFile.Open(args.input, "READ")
if not tfile or tfile.IsZombie():
    sys.exit(f"ERROR: cannot open {args.input}")

tree = tfile.Get(args.tree)
if not tree:
    sys.exit(f"ERROR: tree '{args.tree}' not found in {args.input}")

print(f"Opened {args.input}  ({tree.GetEntries()} entries in '{args.tree}')")

# --- Collect entries that pass selection ---
clusters = []
for entry in tree:
    n_bins = entry.nProfileBins
    if n_bins == 0:
        continue
    ps = entry.profileStart
    if ps <= 0:
        continue
    E = entry.clusterEnergy
    if E < args.energy_min or E > args.energy_max:
        continue

    measured = [entry.measuredProfile[i] for i in range(n_bins)]
    expected = [entry.expectedProfile[i] for i in range(n_bins)]

    clusters.append({
        "energy":   E,
        "nBins":    n_bins,
        "measured": measured,
        "expected": expected,
    })

print(f"{len(clusters)} clusters pass selection cuts")
if not clusters:
    sys.exit("No clusters pass the selection cuts.")

# --- Average profile ---
max_bins = max(cl["nBins"] for cl in clusters)
sum_meas = [0.0] * max_bins
sum_exp = [0.0] * max_bins
counts = [0] * max_bins

for cl in clusters:
    total = sum(cl["measured"]) or 1.0
    for i in range(cl["nBins"]):
        sum_meas[i] += cl["measured"][i] / total
        sum_exp[i] += cl["expected"][i] / total
        counts[i] += 1

t_avg = [i * BIN_WIDTH for i in range(max_bins) if counts[i] > 0]
avg_meas = [sum_meas[i] / counts[i] for i in range(max_bins) if counts[i] > 0]
avg_exp = [sum_exp[i] / counts[i] for i in range(max_bins) if counts[i] > 0]

# --- Plot ---
c = make_transparent_canvas("c_avg", "Average profile", 1000, 650)

y_max = max(max(avg_meas), max(avg_exp)) * 1.25
frame = c.DrawFrame(0, 0, max(t_avg) * 1.05, y_max)
style_frame(frame, "Depth (X_{0})", "Fractional energy deposit / bin")

g_meas = make_graph(t_avg, avg_meas, COLORS["blue"], line_style=1, line_width=3)
g_exp = make_graph(t_avg, avg_exp, COLORS["nominal"], line_style=2, line_width=3)
g_meas.Draw("L SAME")
g_exp.Draw("L SAME")

leg = ROOT.TLegend(0.55, 0.72, 0.90, 0.88)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextColor(TEXT_COLOR)
leg.SetTextSize(0.030)
leg.AddEntry(g_meas, "Mean measured", "l")
leg.AddEntry(g_exp, "Mean expected", "l")
leg.Draw()

label = ROOT.TLatex()
label.SetNDC(True)
label.SetTextColor(TEXT_COLOR)
label.SetTextSize(0.035)
label.SetTextAlign(13)
label.DrawLatex(0.14, 0.97, f"#font[62]{{Average longitudinal profile}} (N = {len(clusters)} clusters)")

os.makedirs(args.output_dir, exist_ok=True)
out_path = os.path.join(args.output_dir, "average_profile.svg")
print(f"Saving SVG to {out_path} ...")
c.SaveAs(out_path)

tfile.Close()
print("Done.")
