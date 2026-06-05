"""
plot_shower_profiles_root.py

Plots measured and expected longitudinal shower energy-deposition curves
from shower_profiles.root using ROOT (PyROOT) as the plotting backend.

Four output canvases saved to a single PDF:
  1. Individual cluster measured vs expected profiles (up to --max-clusters)
  2. Average measured vs average expected profile over all selected clusters
  3. profileStart distribution
  4. profileDiscrepancy distribution

Usage
-----
    python plot_shower_profiles_root.py [--input shower_profiles.root]
                                        [--output shower_profiles.pdf]
                                        [--max-clusters 6]
                                        [--energy-min 0.5]
                                        [--energy-max 1e9]
                                        [--tree ShowerProfiles]
"""

import argparse
import sys
import array
import ctypes

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--input",        default="shower_profiles.root")
parser.add_argument("--tree",         default="ShowerProfiles")
parser.add_argument("--output",       default="shower_profiles.pdf")
parser.add_argument("--max-clusters", type=int,   default=6)
parser.add_argument("--energy-min",   type=float, default=0.5)
parser.add_argument("--energy-max",   type=float, default=1e9)
args = parser.parse_args()

# ── ROOT imports (after arg parsing so --help works without ROOT) ─────────────
import ROOT
ROOT.gROOT.SetBatch(True)   # never open a display
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)
ROOT.gStyle.SetTitleFontSize(0.04)

BIN_WIDTH = 0.5   # radiation lengths per bin (must match PandoraSettings)

# ── Open tree ─────────────────────────────────────────────────────────────────
tfile = ROOT.TFile.Open(args.input, "READ")
if not tfile or tfile.IsZombie():
    sys.exit(f"ERROR: cannot open {args.input}")

tree = tfile.Get(args.tree)
if not tree:
    sys.exit(f"ERROR: tree '{args.tree}' not found in {args.input}")

print(f"Opened {args.input}  ({tree.GetEntries()} entries in '{args.tree}')")

# ── Collect entries that pass selection ───────────────────────────────────────
clusters = []   # list of dicts, one per passing cluster

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
        "energy":       E,
        "profileStart": ps,
        "profileDisc":  entry.profileDiscrepancy,
        "innerLayer":   entry.innerLayer,
        "showerStart":  entry.showerStartLayer,
        "nBins":        n_bins,
        "measured":     measured,
        "expected":     expected,
    })

print(f"{len(clusters)} clusters pass selection cuts")
if not clusters:
    sys.exit("No clusters pass the selection cuts.")

# ── Helper: make a TGraph from lists ─────────────────────────────────────────
def make_graph(x_vals, y_vals, color, line_style=1, line_width=2,
               marker_style=0):
    n = len(x_vals)
    g = ROOT.TGraph(n,
                    array.array('d', x_vals),
                    array.array('d', y_vals))
    g.SetLineColor(color)
    g.SetLineStyle(line_style)
    g.SetLineWidth(line_width)
    g.SetMarkerStyle(marker_style)
    g.SetMarkerSize(5)
    if g.Integral() > 0:
        g.Scale(1/g.Integral())
    return g

# ── Open PDF output ───────────────────────────────────────────────────────────
# ROOT writes multi-page PDFs by opening with filename[, drawing, closing with filename]
pdf_open  = args.output + "["
pdf_close = args.output + "]"

# ── Per-cluster pages: one page per cluster ───────────────────────────────────
c_single = ROOT.TCanvas("c_single", "Per-cluster", 1000, 600)
c_single.SetLeftMargin(0.12)
c_single.SetBottomMargin(0.12)
 
# Open the PDF on the first page
first_page = True
 
print(f"Writing {len(clusters)} individual cluster pages...")
for idx, cl in enumerate(clusters):
    c_single.Clear()
    #print(cl) 
    t = [i * BIN_WIDTH for i in range(cl["nBins"])]
 
    all_y = cl["measured"] + cl["expected"]
    x_max_s = max(t) * 1.10 if t else 1.0
    #y_max_s = max(all_y) * 0.15 if all_y else 1.0
 
 
    g_s_meas = make_graph(t, cl["measured"], ROOT.kBlue+1, line_style=1,
                          line_width=2, marker_style=20)
    g_s_exp  = make_graph(t, cl["expected"], ROOT.kRed+1,  line_style=2, line_width=2)
    #g_s_meas.Draw("L SAME")
    #g_s_exp.Draw("L SAME")

    scale_meas = g_s_meas.Integral()
    scale_exp = g_s_exp.Integral()

    y_max_s = max(ROOT.TMath.MaxElement(g_s_meas.GetN(), g_s_meas.GetY()),
                  ROOT.TMath.MaxElement(g_s_exp.GetN(), g_s_exp.GetY())) * 1.15 # if all_y else 1.0
    print(y_max_s, g_s_meas.GetMaximum(), g_s_exp.GetMaximum())

    frame_s = c_single.DrawFrame(0, 0, x_max_s, y_max_s)
    frame_s.GetXaxis().SetTitle("Depth (X_{0})")
    frame_s.GetYaxis().SetTitle("Fractional energy deposit [% / 100  / bin]")
    frame_s.SetTitle(
        f"Cluster {idx}  E={cl['energy']:.3f} GeV  "
        f"profileStart={cl['profileStart']:.2f} X_{{0}}  "
        f"#chi={cl['profileDisc']:.4f}"
    )

    g_s_meas.Draw("L SAME")
    g_s_exp.Draw("L SAME")

    leg_s = ROOT.TLegend(0.62, 0.72, 0.88, 0.88)
    leg_s.SetBorderSize(0)
    leg_s.SetFillStyle(0)
    leg_s.SetTextSize(0.030)
    leg_s.AddEntry(g_s_meas, "Measured", "l")
    leg_s.AddEntry(g_s_exp,  "Expected", "l")
    leg_s.Draw()
 
    if first_page:
        c_single.Print(pdf_open)   # opens the PDF
        first_page = False
    c_single.Print(args.output)


# ── Canvas 1: individual cluster profiles ─────────────────────────────────────
c1 = ROOT.TCanvas("c1", "Individual cluster profiles", 1000, 600)
c1.SetLeftMargin(0.12)
c1.SetBottomMargin(0.12)

# Pick up to max_clusters spread across the energy range
clusters_sorted = sorted(clusters, key=lambda c: c["energy"])
n_show = min(args.max_clusters, len(clusters_sorted))
step = max(1, len(clusters_sorted) // n_show)
selected = clusters_sorted[::step][:n_show]

# Colour palette
colours = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+2,
           ROOT.kMagenta+1, ROOT.kOrange+1, ROOT.kCyan+2,
           ROOT.kViolet+1, ROOT.kTeal+2]

# Find axis ranges
all_t, all_y = [], []
for cl in selected:
    t = [i * BIN_WIDTH for i in range(cl["nBins"])]
    all_t += t
    all_y += cl["measured"] + cl["expected"]

x_max = max(all_t) * 1.05
#y_max = max(all_y) * 1.15
y_max = 1.0

# Draw frame
frame1 = c1.DrawFrame(0, 0, x_max, y_max)
frame1.GetXaxis().SetTitle("Depth (X_{0})")
frame1.GetYaxis().SetTitle("Fractional energy deposit [% / 100 / bin]")
frame1.SetTitle("Longitudinal shower profiles: measured (solid) vs expected (dashed)")

leg1 = ROOT.TLegend(0.55, 0.55, 0.88, 0.88)
leg1.SetBorderSize(1)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.015)

graphs_keep = []   # keep references alive
for idx, cl in enumerate(selected):
    col = colours[idx % len(colours)]
    t   = [i * BIN_WIDTH for i in range(cl["nBins"])]

    g_meas = make_graph(t, cl["measured"], col, line_style=1, line_width=2)
    g_exp  = make_graph(t, cl["expected"], col, line_style=2, line_width=2)

    g_meas.Draw("L SAME")
    g_exp.Draw("L SAME")
    graphs_keep += [g_meas, g_exp]

    label = f"E={cl['energy']:.2f} GeV  #chi={cl['profileDisc']:.3f}"
    leg1.AddEntry(g_meas, label, "l")

# Add legend entries explaining line styles
g_dummy_solid = make_graph([0],[0], ROOT.kBlack, 1, 2)
g_dummy_dash  = make_graph([0],[0], ROOT.kBlack, 2, 2)
leg1.AddEntry(g_dummy_solid, "Measured",  "l")
leg1.AddEntry(g_dummy_dash,  "Expected",  "l")
graphs_keep += [g_dummy_solid, g_dummy_dash]
leg1.SetBorderSize(0)

leg1.Draw()
c1.Print(pdf_open)   # opens the PDF and writes first page
c1.Print(args.output)

# ── Canvas 2: average measured vs average expected ────────────────────────────
c2 = ROOT.TCanvas("c2", "Average profiles", 1000, 600)
c2.SetLeftMargin(0.12)
c2.SetBottomMargin(0.12)

# Find common max number of bins
max_bins = max(cl["nBins"] for cl in clusters)
sum_meas = [0.0] * max_bins
sum_exp  = [0.0] * max_bins
counts   = [0]   * max_bins

for cl in clusters:
    for i in range(cl["nBins"]):
        sum_meas[i] += cl["measured"][i]
        sum_exp[i]  += cl["expected"][i]
        counts[i]   += 1

t_avg    = [i * BIN_WIDTH for i in range(max_bins) if counts[i] > 0]
avg_meas = [sum_meas[i] / counts[i] for i in range(max_bins) if counts[i] > 0]
avg_exp  = [sum_exp[i]  / counts[i] for i in range(max_bins) if counts[i] > 0]

#y_max2 = max(max(avg_meas), max(avg_exp)) * 1.15
y_max2 = 0.16
frame2 = c2.DrawFrame(0, 0, max(t_avg) * 1.05, y_max2)
frame2.GetXaxis().SetTitle("Depth (X_{0})")
frame2.GetYaxis().SetTitle("Fractional  energy deposit [% / 100 / bin]")
frame2.SetTitle(f"Average longitudinal profile ({len(clusters)} clusters): measured vs expected")

g_avg_meas = make_graph(t_avg, avg_meas, ROOT.kBlue+1, 1, 3)
g_avg_exp  = make_graph(t_avg, avg_exp,  ROOT.kRed+1,  2, 3)
g_avg_meas.Draw("L SAME")
g_avg_exp.Draw("L SAME")

leg2 = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
leg2.SetBorderSize(0)
leg2.SetFillStyle(0)
leg2.SetTextSize(0.02)
leg2.AddEntry(g_avg_meas, "Mean measured", "l")
leg2.AddEntry(g_avg_exp,  "Mean expected", "l")
leg2.Draw()

c2.Print(args.output)

# ── Canvas 3: profileStart distribution ───────────────────────────────────────
c3 = ROOT.TCanvas("c3", "profileStart", 800, 600)
c3.SetLeftMargin(0.12)
c3.SetBottomMargin(0.12)

ps_vals = [cl["profileStart"] for cl in clusters]
h_ps = ROOT.TH1F("h_ps", "Shower profile start;profileStart (X_{0});Clusters / bin",
                 40, 0, max(ps_vals) * 1.1)
h_ps.SetFillColor(ROOT.kBlue-9)
h_ps.SetLineColor(ROOT.kBlue+1)
for v in ps_vals:
    h_ps.Fill(v)
h_ps.Draw("HIST")
c3.Print(args.output)

# ── Canvas 4: profileDiscrepancy distribution ─────────────────────────────────
c4 = ROOT.TCanvas("c4", "profileDiscrepancy", 800, 600)
c4.SetLeftMargin(0.12)
c4.SetBottomMargin(0.12)

pd_vals = [cl["profileDisc"] for cl in clusters]
h_pd = ROOT.TH1F("h_pd", "Profile discrepancy;profileDiscrepancy (normalised);Clusters / bin",
                 40, 0, max(pd_vals) * 1.1)
h_pd.SetFillColor(ROOT.kRed-9)
h_pd.SetLineColor(ROOT.kRed+1)
for v in pd_vals:
    h_pd.Fill(v)
h_pd.Draw("HIST")
c4.Print(args.output)

# ── Close PDF ─────────────────────────────────────────────────────────────────
c4.Print(pdf_close)

tfile.Close()

# ── Summary ───────────────────────────────────────────────────────────────────
import statistics
energies = [cl["energy"] for cl in clusters]
ps_list  = [cl["profileStart"] for cl in clusters]
pd_list  = [cl["profileDisc"]  for cl in clusters]

print(f"\nSummary statistics ({len(clusters)} selected clusters):")
print(f"  Energy range      : {min(energies):.2f} - {max(energies):.2f} GeV")
print(f"  Mean profileStart : {statistics.mean(ps_list):.3f} +/- {statistics.stdev(ps_list):.3f} X0")
print(f"  Mean discrepancy  : {statistics.mean(pd_list):.4f} +/- {statistics.stdev(pd_list):.4f}")
print(f"Saved {args.output}")
