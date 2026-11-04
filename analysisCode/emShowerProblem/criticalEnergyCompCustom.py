"""
criticalEnergyCompCustom.py

Compares Pandora reconstruction of the same photon sample run with two
LongProfileCriticalEnergy settings (8 MeV vs 80 MeV), styled to match
averageProfileCustom.py.

Produces in <output_dir>:
  ec_compare_average_profile.pdf  - mean measured vs mean expected profile for each E_c
  ec_compare_discrepancy.pdf      - profileDiscrepancy distribution for each E_c
and prints, for each E_c:
  - photon ID efficiency (from PandoraPFOs in the reco .slcio)
  - mean profileDiscrepancy (from the ShowerProfiles tree)

Must run inside the mucoll container (needs pyLCIO).

Usage
-----
    python criticalEnergyCompCustom.py [--input-dir ~/mucoll/jobs/output/ec_compare]
                                       [--output-dir customPDFs]
                                       [--energy-min 12.5]
"""

import argparse
import array
import glob
import math
import os
import sys

parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--input-dir",  default=os.path.expanduser("~/mucoll/jobs/output/ec_compare"))
parser.add_argument("--output-dir", default="customPDFs")
parser.add_argument("--tree",       default="ShowerProfiles")
# Half the gun energy: keeps the main photon shower and drops small fragments
parser.add_argument("--energy-min", type=float, default=12.5)
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
from pyLCIO import IOIMPL

BIN_WIDTH = 0.5   # radiation lengths per bin (must match PandoraSettings)

SETTINGS = [
    # tag,    label,                   color key, line style
    ("8MeV",  "E_{c} = 8 MeV",   "blue",    1),
    ("80MeV", "E_{c} = 80 MeV",  "nominal", 2),
]

# --- Color scheme (matches averageProfileCustom.py) ---
COLORS_HEX = {
    "purple":  "#101c4d",
    "nominal": "#e38153",
    "vacuum":  "#d14382",
    "blue":    "#0052c5",
    "gold":    "#e3ab0f",
    "yellow":  "#ffe436",
    "text":    "#000000",
    "axis":    "#000000",
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


def make_canvas(name, title, w, h):
    # Opaque white background so the black text/axes are readable anywhere
    c = ROOT.TCanvas(name, title, w, h)
    c.SetFillColor(ROOT.kWhite)
    c.SetFillStyle(1001)
    c.SetFrameFillColor(ROOT.kWhite)
    c.SetFrameFillStyle(1001)
    c.SetFrameLineColor(AXIS_COLOR)
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
        axis.SetNdivisions(510)
    frame.SetTitle("")


def make_graph(x_vals, y_vals, color, line_style=1, line_width=3):
    n = len(x_vals)
    g = ROOT.TGraph(n, array.array('d', x_vals), array.array('d', y_vals))
    g.SetLineColor(color)
    g.SetLineStyle(line_style)
    g.SetLineWidth(line_width)
    return g


def make_legend(x1, y1, x2, y2):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextColor(TEXT_COLOR)
    leg.SetTextSize(0.030)
    return leg


def draw_label():
    label_top = ROOT.TLatex()
    label_top.SetNDC(True)
    label_top.SetTextColor(TEXT_COLOR)
    label_top.SetTextFont(72)  # Helvetica bold italic
    label_top.SetTextAlign(33)
    label_top.SetTextSize(0.035)
    label_top.DrawLatex(0.93, 0.86, "Muon Collider")

    label_sub = ROOT.TLatex()
    label_sub.SetNDC(True)
    label_sub.SetTextColor(TEXT_COLOR)
    label_sub.SetTextAlign(33)
    label_sub.SetTextSize(0.028)
    label_sub.DrawLatex(0.93, 0.81, "Photon Gun #it{ E}_{#gamma} = 25 GeV")


# --- Photon ID efficiency from the reco .slcio files ---
def photon_id_efficiency(reco_files):
    n_events = n_any_photon = n_leading_photon = 0
    for path in reco_files:
        reader = IOIMPL.LCFactory.getInstance().createLCReader()
        reader.open(path)
        for event in reader:
            n_events += 1
            pfos = list(event.getCollection("PandoraPFOs"))
            if not pfos:
                continue
            if any(p.getType() == 22 for p in pfos):
                n_any_photon += 1
            if max(pfos, key=lambda p: p.getEnergy()).getType() == 22:
                n_leading_photon += 1
        reader.close()
    return n_events, n_any_photon, n_leading_photon


# --- Profiles and discrepancy from the ShowerProfiles tree ---
def read_profiles(profile_files):
    chain = ROOT.TChain(args.tree)
    for path in profile_files:
        chain.Add(path)

    clusters = []
    for entry in chain:
        n_bins = entry.nProfileBins
        # profileStart == 0 is a valid fit (shower starts at the ECAL face); failures are -1
        if n_bins == 0 or entry.profileStart < 0:
            continue
        if entry.clusterEnergy < args.energy_min:
            continue
        clusters.append({
            "nBins":       n_bins,
            "measured":    [entry.measuredProfile[i] for i in range(n_bins)],
            "expected":    [entry.expectedProfile[i] for i in range(n_bins)],
            "discrepancy": entry.profileDiscrepancy,
            "start":       entry.profileStart,
            "pdg":         entry.pdgCode,
        })
    return clusters


def average_profile(clusters):
    max_bins = max(cl["nBins"] for cl in clusters)
    sum_meas = [0.0] * max_bins
    sum_exp = [0.0] * max_bins
    counts = [0] * max_bins

    for cl in clusters:
        total = sum(cl["measured"]) or 1.0
        # Pandora stores the expected profile starting at bin 0, but compares it to the
        # measured profile shifted by the fitted start (expected[i - offset] vs measured[i]),
        # so apply the same shift here before averaging.
        offset = int(round(cl["start"] / BIN_WIDTH))
        for i in range(cl["nBins"]):
            sum_meas[i] += cl["measured"][i] / total
            if i >= offset:
                sum_exp[i] += cl["expected"][i - offset] / total
            counts[i] += 1

    # Plot at bin centres: measured bin i spans [i, i+1) * BIN_WIDTH in X0
    t_avg = [(i + 0.5) * BIN_WIDTH for i in range(max_bins) if counts[i] > 0]
    avg_meas = [sum_meas[i] / counts[i] for i in range(max_bins) if counts[i] > 0]
    avg_exp = [sum_exp[i] / counts[i] for i in range(max_bins) if counts[i] > 0]
    return t_avg, avg_meas, avg_exp


def mean(values):
    return sum(values) / len(values) if values else float("nan")


results = {}
for tag, _, _, _ in SETTINGS:
    reco_files = sorted(glob.glob(os.path.join(args.input_dir, f"reco_ec_compare_{tag}_*.slcio")))
    profile_files = sorted(glob.glob(os.path.join(args.input_dir, f"shower_profiles_ec_compare_{tag}_*.root")))
    if not reco_files or not profile_files:
        sys.exit(f"ERROR: missing reco/profile files for {tag} in {args.input_dir}")

    print(f"[{tag}] {len(reco_files)} reco files, {len(profile_files)} profile files")
    clusters = read_profiles(profile_files)
    if not clusters:
        sys.exit(f"No clusters pass the selection cuts for {tag}.")
    results[tag] = {
        "eff":      photon_id_efficiency(reco_files),
        "clusters": clusters,
        "profile":  average_profile(clusters),
    }

# --- Printout ---
print()
print(f"Clusters with E_EM > {args.energy_min} GeV and a valid longitudinal profile")
print(f"{'':>7} {'events':>7} {'eff(any PFO=gamma)':>20} {'eff(leading PFO=gamma)':>24} "
      f"{'clusters':>9} {'<discrepancy>':>14} {'<discrepancy> PFO=gamma':>24}")
for tag, _, _, _ in SETTINGS:
    n_ev, n_any, n_lead = results[tag]["eff"]
    clusters = results[tag]["clusters"]
    disc_all = [cl["discrepancy"] for cl in clusters if cl["discrepancy"] >= 0]
    disc_gamma = [cl["discrepancy"] for cl in clusters if cl["discrepancy"] >= 0 and cl["pdg"] == 22]
    print(f"{tag:>7} {n_ev:>7d} "
          f"{n_any / n_ev:>14.3f} ± {math.sqrt(n_any / n_ev * (1 - n_any / n_ev) / n_ev):.3f} "
          f"{n_lead / n_ev:>18.3f} ± {math.sqrt(n_lead / n_ev * (1 - n_lead / n_ev) / n_ev):.3f} "
          f"{len(disc_all):>9d} {mean(disc_all):>14.4f} {mean(disc_gamma):>24.4f}")

os.makedirs(args.output_dir, exist_ok=True)

# --- Plot 1: average measured vs expected profile ---
c = make_canvas("c_ec_profile", "Average profile", 1000, 650)
y_max = max(max(max(r["profile"][1]), max(r["profile"][2])) for r in results.values()) * 1.25
t_max = max(max(r["profile"][0]) for r in results.values())
frame = c.DrawFrame(0, 0, t_max * 1.05, y_max)
style_frame(frame, "Depth (X_{0})", "Fractional energy deposit / bin")

# E_c only enters the expected profile, so both runs (same sim hits) give the same
# measured profile. Draw the 8 MeV measured curve as markers so it stays visible on
# top of the 80 MeV line instead of being hidden underneath it.
t8, m8, _ = results["8MeV"]["profile"]
t80, m80, _ = results["80MeV"]["profile"]
n_common = min(len(m8), len(m80))
print(f"Max |measured(8 MeV) - measured(80 MeV)| per bin: "
      f"{max(abs(a - b) for a, b in zip(m8[:n_common], m80[:n_common])):.2e}")

leg = make_legend(0.55, 0.40, 0.92, 0.76)
graphs = []
for tag, label, color, _ in reversed(SETTINGS):
    t_avg, avg_meas, avg_exp = results[tag]["profile"]
    g_meas = make_graph(t_avg, avg_meas, COLORS[color], line_style=1, line_width=3)
    g_exp = make_graph(t_avg, avg_exp, COLORS[color], line_style=2, line_width=3)
    meas_opt = "L SAME"
    if tag == "8MeV":
        g_meas.SetMarkerStyle(24)
        g_meas.SetMarkerSize(0.9)
        g_meas.SetMarkerColor(COLORS[color])
        meas_opt = "P SAME"
    g_meas.Draw(meas_opt)
    g_exp.Draw("L SAME")
    meas_leg_opt = "p" if tag == "8MeV" else "l"
    graphs += [g_meas, g_exp]
    leg.AddEntry(g_meas, f"Mean measured, {label}", meas_leg_opt)
    leg.AddEntry(g_exp, f"Mean expected, {label}", "l")
leg.Draw()
draw_label()

out_path = os.path.join(args.output_dir, "ec_compare_average_profile.pdf")
print(f"Saving PDF to {out_path} ...")
c.SaveAs(out_path)

# --- Plot 2: profileDiscrepancy distribution ---
c2 = make_canvas("c_ec_disc", "Profile discrepancy", 1000, 650)
hists = []
for tag, label, color, line_style in SETTINGS:
    h = ROOT.TH1D(f"h_disc_{tag}", "", 60, 0, 1.5)
    for cl in results[tag]["clusters"]:
        if cl["discrepancy"] >= 0:
            h.Fill(cl["discrepancy"])
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(COLORS[color])
    h.SetLineStyle(line_style)
    h.SetLineWidth(3)
    hists.append((h, label))

frame2 = c2.DrawFrame(0, 0, 1.5, max(h.GetMaximum() for h, _ in hists) * 1.25)
style_frame(frame2, "Longitudinal profile discrepancy", "Fraction of clusters")
leg2 = make_legend(0.16, 0.74, 0.45, 0.88)
for h, label in hists:
    h.Draw("HIST SAME")
    leg2.AddEntry(h, label, "l")
leg2.Draw()
draw_label()

out_path = os.path.join(args.output_dir, "ec_compare_discrepancy.pdf")
print(f"Saving PDF to {out_path} ...")
c2.SaveAs(out_path)
print("Done.")
