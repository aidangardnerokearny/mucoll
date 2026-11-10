"""Average measured profile, nominal vs up_vac, for profileStart > 0 and >= 0.
Styled like averageProfileCustom.py. Usage: python profileStartCutCompCustom.py [out.pdf]  (default: showerProblemPlots/profileStart_cut_nominal_upvac.pdf)"""
import array, glob, os, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

BW, EMIN = 0.5, 0.5
BASE = os.path.expanduser("~/mucoll/jobs/output")
DATASETS = [
    ("Nominal", [f"{BASE}/shower_profiles_nominal.root"], "#e38153"),
    ("Upstream vacuum", sorted(glob.glob(f"{BASE}/up_vac/shower_profiles_up_vac_[0-9]*.root")), "#0052c5"),
]
TEXT = ROOT.kBlack  # all plots use black axes and text


def avg_profile(files, strict):
    ch = ROOT.TChain("ShowerProfiles")
    for f in files:
        ch.Add(f)
    s, c, n = np.zeros(200), np.zeros(200), 0
    for e in ch:
        nb = e.nProfileBins
        if nb == 0 or e.profileStart < 0 or (strict and e.profileStart == 0):
            continue
        if e.clusterEnergy < EMIN:
            continue
        m = np.array([e.measuredProfile[i] for i in range(nb)], dtype=float)
        s[:nb] += m / (m.sum() or 1.0)
        c[:nb] += 1
        n += 1
    k = c > 0
    return (np.arange(200) * BW)[k], (s / np.maximum(c, 1))[k], n


c = ROOT.TCanvas("c", "", 1000, 650)
c.SetFillColorAlpha(ROOT.kWhite, 0); c.SetFillStyle(4000); c.SetFrameFillStyle(4000)
c.SetLeftMargin(0.13); c.SetBottomMargin(0.13); c.SetTopMargin(0.10); c.SetRightMargin(0.05)
c.SetTickx(1); c.SetTicky(1)

curves = []
for name, files, hexcol in DATASETS:
    col = ROOT.TColor.GetColor(hexcol)
    for strict, style, lab in ((False, 1, "#geq 0"), (True, 2, "> 0")):
        t, y, n = avg_profile(files, strict)
        curves.append((f"{name}, profileStart {lab} (N = {n})", t, y, col, style))

ymax = max(y.max() for _, _, y, _, _ in curves) * 1.25
fr = c.DrawFrame(0, 0, 30, ymax)
fr.SetTitle("")
fr.GetXaxis().SetTitle("Depth (X_{0})")
fr.GetYaxis().SetTitle("Fractional energy deposit / bin")
for ax in (fr.GetXaxis(), fr.GetYaxis()):
    ax.SetTitleColor(TEXT); ax.SetLabelColor(TEXT); ax.SetAxisColor(TEXT)
    ax.SetTickLength(0.02); ax.SetNdivisions(510)

leg = ROOT.TLegend(0.45, 0.45, 0.92, 0.72)
leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextColor(TEXT); leg.SetTextSize(0.028)
graphs = []
for lab, t, y, col, style in curves:
    g = ROOT.TGraph(len(t), array.array("d", t), array.array("d", y))
    g.SetLineColor(col); g.SetLineStyle(style); g.SetLineWidth(3); g.Draw("L SAME")
    leg.AddEntry(g, lab, "l"); graphs.append(g)
leg.Draw()

lt = ROOT.TLatex(); lt.SetNDC(True); lt.SetTextColor(TEXT); lt.SetTextAlign(33)
lt.SetTextFont(72); lt.SetTextSize(0.035); lt.DrawLatex(0.93, 0.86, "Muon Collider")
lt.SetTextFont(42); lt.SetTextSize(0.028); lt.DrawLatex(0.93, 0.81, "Photon Gun #it{ E}_{#gamma} = 25 GeV")

out = sys.argv[1] if len(sys.argv) > 1 else "showerProblemPlots/profileStart_cut_nominal_upvac.pdf"
c.SaveAs(out)
for lab, t, y, *_ in curves:
    print(f"{lab:55s} peak {t[np.argmax(y)]:.1f} X0, mean {(t*y).sum()/y.sum():.2f} X0")
