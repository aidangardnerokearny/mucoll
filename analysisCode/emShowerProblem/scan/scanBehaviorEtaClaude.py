#!/usr/bin/env python3
"""
Average per-cluster shower observables as a function of polar angle theta,
for the nominal geometry vs. the vacuum-solenoid geometry.

One page per observable: top pad overlays <var>(theta) for both geometries,
bottom pad shows the vacuum - nominal difference.

Errors are the standard error on the mean in each theta bin (RMS/sqrt(N)),
which is what TProfile stores with the default error option. The difference
pad is computed *per matched entry* when the two trees are entry-aligned, so
the error on the difference automatically accounts for the correlation
between the two geometries (same event, same shower); it falls back to
adding the two errors in quadrature if the trees are not aligned.
"""

import math
import argparse
import sys
import array

# --- Arguments ---
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--nominal",       default="shower_profiles.root")
parser.add_argument("--vacuum",        default="shower_profiles_vac.root")
parser.add_argument("--tree",          default="ShowerProfiles")
parser.add_argument("--output",        default="scan_theta.pdf")
parser.add_argument("--energy-min",    type=float, default=0.5)
parser.add_argument("--energy-max",    type=float, default=1e9)
parser.add_argument("--theta-branch",  default=None,
                    help="Branch holding theta. Default: auto-detect.")
parser.add_argument("--theta-units",   default="auto", choices=["auto", "rad", "deg"])
parser.add_argument("--theta-min",     type=float, default=0.0,   help="degrees")
parser.add_argument("--theta-max",     type=float, default=180.0, help="degrees")
parser.add_argument("--theta-bins",    type=int,   default=36)
parser.add_argument("--min-entries",   type=int,   default=5,
                    help="Skip theta bins with fewer clusters than this.")
parser.add_argument("--fold",          action="store_true",
                    help="Fold theta -> 180-theta for theta > 90 (forward/backward symmetry).")
args = parser.parse_args()

# --- ROOT imports ---
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)
ROOT.gStyle.SetTitleFontSize(0.04)

COL_NOM = ROOT.kBlue + 1
COL_VAC = ROOT.kRed + 1
COL_DIF = ROOT.kBlack

# Observables that are meaningful for every cluster (no profile required).
BASIC_VARS = {"clusterEnergy", "nCaloHits"}

# Axis labels / units for the branches we expect. Anything not listed here is
# still plotted, just with the raw branch name as the label.
LABELS = {
    "clusterEnergy":       ("Cluster energy", "GeV"),
    "nCaloHits":           ("Number of calo hits", ""),
    "nProfileBins":        ("Number of profile bins", ""),
    "profileStart":        ("Profile start", "X_{0}"),
    "profileDiscrepancy":  ("Profile discrepancy", ""),
    "innerLayer":          ("Inner pseudo-layer", ""),
    "showerStartLayer":    ("Shower start pseudo-layer", ""),
    "peakRms":             ("Peak RMS", ""),
    "PeakRms":             ("Peak RMS", ""),
    "rmsRatio":            ("RMS ratio", ""),
    "RmsRatio":            ("RMS ratio", ""),
}

# Candidates used to work out theta if --theta-branch is not given.
THETA_CANDIDATES = ["clusterTheta", "theta", "clusterDirTheta", "mcTheta"]
VECTOR_CANDIDATES = [
    ("clusterDirX",       "clusterDirY",       "clusterDirZ"),
    ("clusterDirectionX", "clusterDirectionY", "clusterDirectionZ"),
    ("clusterCoGX",       "clusterCoGY",       "clusterCoGZ"),
    ("centroidX",         "centroidY",         "centroidZ"),
    ("mcMomentumX",       "mcMomentumY",       "mcMomentumZ"),
]


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

trees = [tree_nom, tree_vac]
tags = ["Nominal", "Vacuum"]

print(f"Opened {args.nominal}: {tree_nom.GetEntries()} entries in {args.tree}")
print(f"Opened {args.vacuum}: {tree_vac.GetEntries()} entries in {args.tree}")

paired = (tree_nom.GetEntries() == tree_vac.GetEntries())
if not paired:
    print("WARNING: the two trees do not have the same number of entries; "
          "the difference pad will use uncorrelated (quadrature) errors.")


# --- Work out which branches are scalars, and where theta comes from ---
def branch_names(tree):
    return [b.GetName() for b in tree.GetListOfBranches()]


def scalar_branches(tree):
    """Branch names holding a single number per entry (no variable-length arrays)."""
    out = []
    for br in tree.GetListOfBranches():
        leaves = br.GetListOfLeaves()
        if leaves.GetEntries() != 1:
            continue
        leaf = leaves.At(0)
        if leaf.GetLeafCount():          # variable-length array, e.g. measuredProfile
            continue
        if leaf.GetLen() != 1:           # fixed-length array
            continue
        out.append(br.GetName())
    return out


def resolve_theta(tree):
    """Return (mode, spec): ('branch', (name,)) or ('vector', (x, y, z))."""
    names = set(branch_names(tree))

    if args.theta_branch:
        if args.theta_branch not in names:
            sys.exit(f"ERROR: branch '{args.theta_branch}' not in tree '{args.tree}'")
        return "branch", (args.theta_branch,)

    for name in THETA_CANDIDATES:
        if name in names:
            return "branch", (name,)

    for trip in VECTOR_CANDIDATES:
        if all(c in names for c in trip):
            return "vector", trip

    sys.exit("ERROR: could not find a theta branch or an (x, y, z) triplet.\n"
             f"       Branches available: {sorted(names)}\n"
             "       Pass one explicitly with --theta-branch.")


theta_mode, theta_spec = resolve_theta(tree_nom)
if theta_mode == "branch":
    print(f"Using theta from branch '{theta_spec[0]}'")
else:
    print(f"Computing theta from ({', '.join(theta_spec)})")

scalars_nom = scalar_branches(tree_nom)
scalars_vac = set(scalar_branches(tree_vac))
plot_vars = [v for v in scalars_nom
             if v in scalars_vac and v not in theta_spec]

if not plot_vars:
    sys.exit("ERROR: no scalar branches in common between the two trees")

print(f"Plotting {len(plot_vars)} observables: {', '.join(plot_vars)}")


def get_theta_raw(event):
    if theta_mode == "branch":
        return getattr(event, theta_spec[0])
    x = getattr(event, theta_spec[0])
    y = getattr(event, theta_spec[1])
    z = getattr(event, theta_spec[2])
    r = math.sqrt(x * x + y * y)
    if r == 0.0 and z == 0.0:
        return None
    return math.atan2(r, z)


# --- Collect entries ---
def collect(tree):
    """One record per entry, keeping the selection flags rather than dropping rows,
    so that the nominal and vacuum lists stay index-aligned."""
    records = []
    for event in tree:
        theta_raw = get_theta_raw(event)
        if theta_raw is None:
            records.append(None)
            continue

        vals = {v: float(getattr(event, v)) for v in plot_vars}

        E = vals.get("clusterEnergy", 0.0)
        pass_basic = (args.energy_min <= E <= args.energy_max)

        n_bins = int(vals.get("nProfileBins", 0))
        ps = vals.get("profileStart", 0.0)
        pass_profile = pass_basic and n_bins > 0 and ps > 0.0

        records.append({
            "thetaRaw":    theta_raw,
            "vals":        vals,
            "passBasic":   pass_basic,
            "passProfile": pass_profile,
        })
    return records


records_nom = collect(tree_nom)
records_vac = collect(tree_vac)
both_records = [records_nom, records_vac]

# --- Radians or degrees? ---
raw_all = [r["thetaRaw"] for r in records_nom + records_vac if r]
if not raw_all:
    sys.exit("ERROR: no usable theta values")

if args.theta_units == "auto":
    is_rad = max(abs(v) for v in raw_all) <= 1.05 * math.pi
else:
    is_rad = (args.theta_units == "rad")
print(f"Interpreting theta as {'radians' if is_rad else 'degrees'}")

for records in both_records:
    for r in records:
        if not r:
            continue
        theta = math.degrees(r["thetaRaw"]) if is_rad else r["thetaRaw"]
        if args.fold and theta > 90.0:
            theta = 180.0 - theta
        r["theta"] = theta

for records, tag in zip(both_records, tags):
    n_basic = sum(1 for r in records if r and r["passBasic"])
    n_prof = sum(1 for r in records if r and r["passProfile"])
    print(f"{tag}: {n_basic} clusters pass the energy cut, "
          f"{n_prof} also have a valid profile")
    if n_basic == 0:
        sys.exit(f"No {tag} clusters pass the selection cuts")


# --- Helpers ---
def label_for(var):
    name, unit = LABELS.get(var, (var, ""))
    return f"#LT{name}#GT [{unit}]" if unit else f"#LT{name}#GT"


def selected(records, need_profile):
    key = "passProfile" if need_profile else "passBasic"
    for r in records:
        if r and r[key]:
            yield r


def make_profile(name, records, var, need_profile):
    """TProfile of <var> vs theta. Default error option => RMS/sqrt(N) per bin."""
    prof = ROOT.TProfile(name, "", args.theta_bins, args.theta_min, args.theta_max)
    prof.SetDirectory(0)
    for r in selected(records, need_profile):
        prof.Fill(r["theta"], r["vals"][var])
    return prof


def make_diff_profile(name, var, need_profile):
    """Per-entry (vacuum - nominal) vs theta, using entries that pass in both trees."""
    prof = ROOT.TProfile(name, "", args.theta_bins, args.theta_min, args.theta_max)
    prof.SetDirectory(0)
    key = "passProfile" if need_profile else "passBasic"
    for r_nom, r_vac in zip(records_nom, records_vac):
        if not r_nom or not r_vac:
            continue
        if not (r_nom[key] and r_vac[key]):
            continue
        prof.Fill(r_nom["theta"], r_vac["vals"][var] - r_nom["vals"][var])
    return prof


def make_graph(prof, color, line_style=1, line_width=2, marker_style=20):
    """TGraphErrors from the filled bins of a TProfile."""
    x_vals, y_vals, ex_vals, ey_vals = [], [], [], []
    for b in range(1, prof.GetNbinsX() + 1):
        if prof.GetBinEntries(b) < args.min_entries:
            continue
        x_vals.append(prof.GetXaxis().GetBinCenter(b))
        y_vals.append(prof.GetBinContent(b))
        ex_vals.append(0.5 * prof.GetXaxis().GetBinWidth(b))
        ey_vals.append(prof.GetBinError(b))

    g = ROOT.TGraphErrors(len(x_vals),
                          array.array('d', x_vals),
                          array.array('d', y_vals),
                          array.array('d', ex_vals),
                          array.array('d', ey_vals))
    g.SetLineColor(color)
    g.SetLineStyle(line_style)
    g.SetLineWidth(line_width)
    g.SetMarkerColor(color)
    g.SetMarkerStyle(marker_style)
    g.SetMarkerSize(0.9)
    return g


def quadrature_diff(g_nom, g_vac):
    """Fallback difference graph when the trees are not entry-aligned."""
    y_nom = {round(g_nom.GetPointX(i), 6): i for i in range(g_nom.GetN())}
    x_vals, y_vals, ex_vals, ey_vals = [], [], [], []
    for j in range(g_vac.GetN()):
        x = round(g_vac.GetPointX(j), 6)
        if x not in y_nom:
            continue
        i = y_nom[x]
        x_vals.append(x)
        y_vals.append(g_vac.GetPointY(j) - g_nom.GetPointY(i))
        ex_vals.append(g_vac.GetErrorX(j))
        ey_vals.append(math.hypot(g_vac.GetErrorY(j), g_nom.GetErrorY(i)))

    g = ROOT.TGraphErrors(len(x_vals),
                          array.array('d', x_vals),
                          array.array('d', y_vals),
                          array.array('d', ex_vals),
                          array.array('d', ey_vals))
    g.SetLineColor(COL_DIF)
    g.SetLineWidth(2)
    g.SetMarkerColor(COL_DIF)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.9)
    return g


def y_range(graphs, pad=0.15):
    lo, hi = None, None
    for g in graphs:
        for i in range(g.GetN()):
            y, e = g.GetPointY(i), g.GetErrorY(i)
            lo = y - e if lo is None else min(lo, y - e)
            hi = y + e if hi is None else max(hi, y + e)
    if lo is None:
        return 0.0, 1.0
    if hi == lo:
        return lo - 1.0, hi + 1.0
    span = hi - lo
    return lo - pad * span, hi + pad * span


# --- Open PDF output ---
pdf_open = args.output + "["
pdf_close = args.output + "]"

c = ROOT.TCanvas("c_theta", "vs theta", 1000, 700)
c.Print(pdf_open)

theta_title = "#theta [deg]" if not args.fold else "#theta (folded) [deg]"
keep = []   # keep ROOT objects alive until the PDF is closed

# --- Page 1: occupancy, so the reader can see where the statistics are ---
c.Clear()
c.SetLeftMargin(0.12)
c.SetBottomMargin(0.12)

h_nom = ROOT.TH1D("h_theta_nom", f"Cluster occupancy;{theta_title};Clusters",
                  args.theta_bins, args.theta_min, args.theta_max)
h_vac = ROOT.TH1D("h_theta_vac", "", args.theta_bins, args.theta_min, args.theta_max)
h_nom.SetDirectory(0)
h_vac.SetDirectory(0)
for r in selected(records_nom, True):
    h_nom.Fill(r["theta"])
for r in selected(records_vac, True):
    h_vac.Fill(r["theta"])
h_nom.SetLineColor(COL_NOM)
h_nom.SetLineWidth(2)
h_vac.SetLineColor(COL_VAC)
h_vac.SetLineWidth(2)
h_vac.SetLineStyle(2)
h_nom.SetMaximum(1.25 * max(h_nom.GetMaximum(), h_vac.GetMaximum(), 1.0))
h_nom.Draw("HIST")
h_vac.Draw("HIST SAME")

leg0 = ROOT.TLegend(0.62, 0.75, 0.88, 0.88)
leg0.SetBorderSize(0)
leg0.SetFillStyle(0)
leg0.AddEntry(h_nom, "Nominal", "l")
leg0.AddEntry(h_vac, "Vacuum solenoid", "l")
leg0.Draw()
keep += [h_nom, h_vac, leg0]
c.Print(args.output)

# --- One page per observable ---
for var in plot_vars:
    need_profile = var not in BASIC_VARS

    p_nom = make_profile(f"p_nom_{var}", records_nom, var, need_profile)
    p_vac = make_profile(f"p_vac_{var}", records_vac, var, need_profile)
    g_nom = make_graph(p_nom, COL_NOM, marker_style=20)
    g_vac = make_graph(p_vac, COL_VAC, line_style=2, marker_style=24)

    if g_nom.GetN() == 0 and g_vac.GetN() == 0:
        print(f"Skipping {var}: no theta bin has >= {args.min_entries} clusters")
        continue

    if paired:
        g_dif = make_graph(make_diff_profile(f"p_dif_{var}", var, need_profile),
                           COL_DIF, marker_style=20)
        dif_note = "per-event difference"
    else:
        g_dif = quadrature_diff(g_nom, g_vac)
        dif_note = "quadrature errors"

    c.Clear()
    pad_top = ROOT.TPad(f"top_{var}", "", 0.0, 0.30, 1.0, 1.0)
    pad_bot = ROOT.TPad(f"bot_{var}", "", 0.0, 0.00, 1.0, 0.30)
    pad_top.SetBottomMargin(0.02)
    pad_top.SetLeftMargin(0.12)
    pad_top.SetTicks(1, 1)
    pad_bot.SetTopMargin(0.04)
    pad_bot.SetBottomMargin(0.32)
    pad_bot.SetLeftMargin(0.12)
    pad_bot.SetTicks(1, 1)
    pad_top.Draw()
    pad_bot.Draw()

    # top pad: the two means
    pad_top.cd()
    sample = "profile-selected clusters" if need_profile else "all clusters"
    lo, hi = y_range([g_nom, g_vac])
    frame_top = pad_top.DrawFrame(args.theta_min, lo, args.theta_max, hi)
    frame_top.SetTitle(f"{LABELS.get(var, (var, ''))[0]} vs #theta  ({sample})")
    frame_top.GetYaxis().SetTitle(label_for(var))
    frame_top.GetYaxis().SetTitleSize(0.055)
    frame_top.GetYaxis().SetTitleOffset(1.05)
    frame_top.GetYaxis().SetLabelSize(0.045)
    frame_top.GetXaxis().SetLabelSize(0.0)
    g_nom.Draw("P SAME")
    g_vac.Draw("P SAME")

    n_nom = sum(1 for _ in selected(records_nom, need_profile))
    n_vac = sum(1 for _ in selected(records_vac, need_profile))
    leg = ROOT.TLegend(0.58, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(g_nom, f"Nominal (N = {n_nom})", "lp")
    leg.AddEntry(g_vac, f"Vacuum solenoid (N = {n_vac})", "lp")
    leg.Draw()

    # bottom pad: vacuum - nominal
    pad_bot.cd()
    dlo, dhi = y_range([g_dif], pad=0.25)
    dlo, dhi = min(dlo, 0.0), max(dhi, 0.0)
    frame_bot = pad_bot.DrawFrame(args.theta_min, dlo, args.theta_max, dhi)
    frame_bot.GetXaxis().SetTitle(theta_title)
    frame_bot.GetYaxis().SetTitle("Vac #minus Nom")
    frame_bot.GetXaxis().SetTitleSize(0.13)
    frame_bot.GetXaxis().SetLabelSize(0.11)
    frame_bot.GetXaxis().SetTitleOffset(1.05)
    frame_bot.GetYaxis().SetTitleSize(0.11)
    frame_bot.GetYaxis().SetTitleOffset(0.45)
    frame_bot.GetYaxis().SetLabelSize(0.10)
    frame_bot.GetYaxis().SetNdivisions(505)

    zero = ROOT.TLine(args.theta_min, 0.0, args.theta_max, 0.0)
    zero.SetLineStyle(2)
    zero.SetLineColor(ROOT.kGray + 2)
    zero.Draw()
    g_dif.Draw("P SAME")

    note = ROOT.TLatex()
    note.SetNDC()
    note.SetTextSize(0.10)
    note.SetTextColor(ROOT.kGray + 2)
    note.DrawLatex(0.60, 0.90, dif_note)

    keep += [pad_top, pad_bot, g_nom, g_vac, g_dif, leg, zero, note,
             p_nom, p_vac]
    c.Print(args.output)
    print(f"  wrote page for {var}")

c.Print(pdf_close)
print(f"Wrote {args.output}")
