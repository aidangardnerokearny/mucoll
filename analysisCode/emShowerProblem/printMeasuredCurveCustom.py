"""
printMeasuredCurveCustom.py

PyROOT script, styled to match averageProfileCustom.py / singleClusterProfileCustom.py,
that plots the measured longitudinal shower profile for one event/cluster
(the --cluster-index'th cluster passing the selection cuts).

Produces a single SVG: <output_dir>/measured_curve.svg
  - measured longitudinal shower profile for the selected cluster.

Usage
-----
    python printMeasuredCurveCustom.py [--input shower_profiles.root]
                                       [--tree ShowerProfiles]
                                       [--output-dir customSVGs]
                                       [--energy-min 0.5]
                                       [--energy-max 1e9]
                                       [--cluster-index 0]
"""

import argparse
import os
import sys
import array

parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--input",         default="shower_profiles.root")
parser.add_argument("--tree",          default="ShowerProfiles")
parser.add_argument("--output-dir",    default="customSVGs")
parser.add_argument("--energy-min",    type=float, default=0.5)
parser.add_argument("--energy-max",    type=float, default=1e9)
parser.add_argument("--cluster-index", type=int,   default=0,
                     help="index (0-based) among clusters passing the selection cuts")
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

    clusters.append({
        "energy":   E,
        "nBins":    n_bins,
        "measured": measured,
    })

print(f"{len(clusters)} clusters pass selection cuts")
if not clusters:
    sys.exit("No clusters pass the selection cuts.")
if args.cluster_index < 0 or args.cluster_index >= len(clusters):
    sys.exit(f"ERROR: --cluster-index {args.cluster_index} out of range "
              f"(0..{len(clusters) - 1})")

cluster = clusters[args.cluster_index]
energy = cluster["energy"]
n_bins = cluster["nBins"]
measured = cluster["measured"]

t_vals = [i * BIN_WIDTH for i in range(n_bins)]

print(f"Cluster {args.cluster_index} (E = {energy:.2f} GeV, {n_bins} bins)")
print(f"{'Depth (X0)':>12}  {'Measured':>14}")
for t, m in zip(t_vals, measured):
    print(f"{t:12.2f}  {m:14.6f}")

# --- Plot ---
c = make_transparent_canvas("c_measured", "Measured profile", 1000, 650)

y_max = max(measured) * 1.25
frame = c.DrawFrame(0, 0, max(t_vals) * 1.05, y_max)
style_frame(frame, "Depth (X_{0})", "Energy deposit / bin")

g_meas = make_graph(t_vals, measured, COLORS["blue"], line_style=1, line_width=3)
g_meas.Draw("L SAME")

leg = ROOT.TLegend(0.55, 0.12, 0.90, 0.24)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextColor(TEXT_COLOR)
leg.SetTextSize(0.030)
leg.AddEntry(g_meas, "Measured", "l")
leg.Draw()

label_top = ROOT.TLatex()
label_top.SetNDC(True)
label_top.SetTextColor(TEXT_COLOR)
label_top.SetTextFont(72)  # Helvetica bold italic
label_top.SetTextAlign(33)  # right, top
label_top.SetTextSize(0.035)
label_top.DrawLatex(0.93, 0.86, "Muon Collider")

label_sub = ROOT.TLatex()
label_sub.SetNDC(True)
label_sub.SetTextColor(TEXT_COLOR)
label_sub.SetTextAlign(33)
label_sub.SetTextSize(0.028)
label_sub.DrawLatex(0.93, 0.81, f"Photon Gun #it{{ E}}_{{#gamma}} = {energy:.0f} GeV")

os.makedirs(args.output_dir, exist_ok=True)
out_path = os.path.join(args.output_dir, "measured_curve.svg")
print(f"Saving SVG to {out_path} ...")
c.SaveAs(out_path)

tfile.Close()
print("Done.")
