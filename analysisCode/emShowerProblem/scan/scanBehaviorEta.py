import math
import argparse
import sys
import array
import ctypes

# --- Aguments ---
parser = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--nominal",       default="shower_pofiles.root")
parser.add_argument("--vacuum",        default="shower_profiles_vac.root")
parser.add_argument("--tree",          default="ShowerProfiles")
parser.add_argument("--output",        default="scan_showers.pdf")
parser.add_argument("--max-clusters",  type=int,   default=6)
parser.add_argument("--energy-min",    type=float, default=0.5)
parser.add_argument("--energy-max",    type=float, default=1e9)
args = parser.parse_args()

# --- ROOT imports ---
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)
ROOT.gStyle.SetTitleFontSize(0.04)

BIN_WIDTH = 0.5


# --- Read trees ---
tfile_nom = ROOT.TFile.Open(args.nominal, "READ")
tfile_vac = ROOT.TFile.Open(args.vacuum,  "READ")

if not tfile_nom or tfile_nom.IsZombie():
    sys.exit(f"ERROR: Cannot open {args.nominal}")

if not tfile_vac or tfile_vac.IsZombie():
    sys.exit(f"ERROR: Cannot open {args.vacuum}")

tree_nom = tfile_nom.Get(args.tree)
if not tree_nom:
    sys.exit(f"ERROR: Tree '{args.tree}' not found in {args.nominal}")

tree_vac = tfile_vac.Get(args.tree)
if not tree_vac:
    sys.exit(f"ERROR: Tree '{args.tree}' not found in {args.vacuum}")

trees = [tree_nom, tree_v]

print(f"Opened {args.nominal}: {tree_nom.GetEntries()} entries in {args.tree}")
print(f"Opened {args.vacuum}: {tree_vac.GetEntries()} entries in {args.tree}")

if tree_nom.GetEntries() != tree_vac.GetEntries():
    sys.exit(f"""ERROR: {args.nominal} tree {args.tree} does not have the same
             number of entries as {args.vacuum} tree {args.tree}""")

# --- Collect entries that pass selection --- 
clusters_nom = []
all_clusters_nom = []
clusters_vac = []
all_clusters_vac = []
both_clusters = [clusters_nom, clusters_vac]
both_all_clusters = [all_clusters_nom, all_clusters_vac]

for i, tree in enumerate(trees):
    clusters = []
    all_clusters = []

    for event in tree:
        E = entry.clusterEnergy
        n_hits = event.nCaloHits

        all_clusters.append({"energy": E, "nCaloHits": n_hits})

        n_bins = event.nProfileBins
        if n_bins == 0:
            continue
        ps = event.profileStart
        if ps <= 0:
            continue
        E = event.clusterEnergy
        if E < args.energy_min or E > args.energy_max:
            continue

        measured = [event.measuredProfile[i] for i in range(n_bins)]
        expected = [event.expectedProfile[i] for i in range(n_bins)]

        clusters.append({
            "energy":         E,
            "profileStart":   ps,
            "profileDisc":    event.profileDiscrepancy,
            "innerLayer":     event.innerLayer,
            "showerStart":    event.showerStartLayer,
            "nBins":          n_bins,
            "measured":       measured,
            "expected":       expected,
        })

    both_clusters[i] = clusters.copy()
    both_all_clusters[i] = all_clusters.copy()
    print(f"{len(clusters)} clusters pass selection cuts")
    if not clusters:
        sys.exit("No clusters pass the selection cuts")


# --- Helper ---
def make_graph(tx_vals, y_vals, color, line_style=1, line_width=2,
               marker_style=0):
    n = len(x_vals)
    g = ROOT.TGraph(n,
                   array.array('d', x_vals),
                   array.array('d', y_vals))
    g.SetLineColor(color)
    g.SetLineStyle(line_style)
    g.SetLineWidth(line_width)
    g.SetMarkerStyle(marker_style)
    g.SetmarkerSize(5)
    return g


# --- Open PDF output ---
pdf_open = args.output + "["
pdf_close = args.output + "]"

c_single = ROOT.TCanvas("c_single", "Per-Cluster", 1000, 600)
c_single.SetLeftMargin(0.12)
c_single.SetBottomMargin(0.12)

first_page = True
print(f"Writing {len(both_clusters[0])}")

for idx, cl in enumertae(both_clusters[0]):
    pc_single.Clear()
    t = [i * BIN_WIDTH for i in range(cl["nBins"])]

    all_y = cl["Nominal"] + cl["Vacuum"]
    x_max_s = max(t) *1.10 if t else 1.0
    g_s_nom = make_graph(t, cl["nominal"], ROOT.kBlue+1)
    g_s_vac = make_graph(t, cl["nominal"], ROOT.kRed+1, line_style=2)

    scale_nom = g_s_nom.Integral()
    scale_vac = g_s_vac.Integral()
