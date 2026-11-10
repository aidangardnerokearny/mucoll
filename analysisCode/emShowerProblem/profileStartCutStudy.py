"""Compare averaged longitudinal profile for profileStart > 0 vs >= 0, per dataset.
Usage: python profileStartCutStudy.py [output_dir]  (default: showerProblemPlots)"""
import glob, os, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

BW = 0.5
EMIN = 0.5
BASE = os.path.expanduser("~/mucoll/jobs/output")
OUT = sys.argv[1] if len(sys.argv) > 1 else "showerProblemPlots"

DATASETS = {
    "nominal":       [f"{BASE}/shower_profiles_nominal.root"],
    "vacuum":        [f"{BASE}/shower_profiles_vacuum.root"],
    "ec_80MeV":      sorted(glob.glob(f"{BASE}/ec_compare/shower_profiles_ec_compare_80MeV_*.root")),
    "ec_8MeV":       sorted(glob.glob(f"{BASE}/ec_compare/shower_profiles_ec_compare_8MeV_*.root")),
    "up_vac":        sorted(glob.glob(f"{BASE}/up_vac/shower_profiles_up_vac_[0-9]*.root")),
}


def load(files):
    ch = ROOT.TChain("ShowerProfiles")
    for f in files:
        ch.Add(f)
    rows = []
    n_total = n_fail = n_lowE = 0
    for e in ch:
        n_total += 1
        n = e.nProfileBins
        if n == 0 or e.profileStart < 0:
            n_fail += 1
            continue
        if e.clusterEnergy < EMIN:
            n_lowE += 1
            continue
        m = np.array([e.measuredProfile[i] for i in range(n)], dtype=float)
        x = np.array([e.expectedProfile[i] for i in range(n)], dtype=float)
        rows.append(dict(ps=e.profileStart, E=e.clusterEnergy, disc=e.profileDiscrepancy,
                         meas=m, exp=x, inner=e.innerLayer - e.ecalInnerPseudoLayer,
                         pdg=e.pdgCode))
    return rows, n_total, n_fail, n_lowE


def avg(rows, key="meas"):
    if not rows:
        return np.array([]), np.array([])
    nb = max(len(r[key]) for r in rows)
    s = np.zeros(nb); c = np.zeros(nb)
    for r in rows:
        tot = r["meas"].sum() or 1.0
        s[: len(r[key])] += r[key] / tot
        c[: len(r[key])] += 1
    k = c > 0
    return (np.arange(nb) * BW)[k], (s / np.maximum(c, 1))[k]


def shape_stats(t, y):
    y = y / y.sum()
    mean = (t * y).sum()
    return t[np.argmax(y)], mean


summary = []
canvases = []
pdf = os.path.join(OUT, "profileStart_cut_study.pdf")
c = ROOT.TCanvas("c", "", 1100, 500)
c.Print(pdf + "[")
for name, files in DATASETS.items():
    if not files:
        continue
    rows, nt, nf, nl = load(files)
    gt = [r for r in rows if r["ps"] > 0]
    eq = [r for r in rows if r["ps"] == 0]
    ge = rows
    res = dict(name=name, total=nt, fail=nf, lowE=nl, n_ge=len(ge), n_gt=len(gt), n_eq=len(eq))
    curves = {}
    for lab, sel in (("ps>0", gt), ("ps>=0", ge), ("ps==0", eq)):
        t, y = avg(sel)
        curves[lab] = (t, y)
        if len(t):
            pk, mn = shape_stats(t, y)
            res[lab] = dict(peak=pk, mean=mn,
                            E=np.mean([r["E"] for r in sel]),
                            disc=np.mean([r["disc"] for r in sel]),
                            inner=np.mean([r["inner"] for r in sel]),
                            first_bin=y[0],
                            photon=np.mean([r["pdg"] == 22 for r in sel]))
    # ps distribution
    res["ps_hist"] = np.histogram([r["ps"] for r in rows], bins=np.arange(-0.25, 10.25, 0.5))[0]
    summary.append(res)

    # plot: left = overlay, right = ratio (>=0)/(>0)
    c.Clear(); c.Divide(2, 1)
    c.cd(1)
    ymax = max(curves[k][1].max() for k in curves if len(curves[k][0])) * 1.2
    xmax = max(curves[k][0].max() for k in curves if len(curves[k][0]))
    fr = ROOT.gPad.DrawFrame(0, 0, min(xmax, 40), ymax)
    fr.SetTitle(f"{name};Depth (X_{{0}});Fractional energy / bin")
    leg = ROOT.TLegend(0.45, 0.65, 0.88, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0)
    gs = []
    for lab, col, sty in (("ps>0", ROOT.kBlue + 1, 1), ("ps>=0", ROOT.kRed + 1, 2), ("ps==0", ROOT.kGray + 1, 3)):
        t, y = curves[lab]
        if not len(t):
            continue
        g = ROOT.TGraph(len(t), t.astype("d"), y.astype("d"))
        g.SetLineColor(col); g.SetLineStyle(sty); g.SetLineWidth(3); g.Draw("L SAME")
        n = {"ps>0": res["n_gt"], "ps>=0": res["n_ge"], "ps==0": res["n_eq"]}[lab]
        leg.AddEntry(g, f"{lab.replace('ps', 'profileStart ')} (N={n})", "l"); gs.append(g)
    leg.Draw()
    c.cd(2)
    t1, y1 = curves["ps>0"]; t2, y2 = curves["ps>=0"]
    n = min(len(t1), len(t2))
    r = np.where(y1[:n] > 0, y2[:n] / np.where(y1[:n] > 0, y1[:n], 1), 1.0)
    fr2 = ROOT.gPad.DrawFrame(0, 0.5, min(xmax, 40), 1.5)
    fr2.SetTitle(f"{name};Depth (X_{{0}});(>= 0) / (> 0)")
    gr = ROOT.TGraph(n, t1[:n].astype("d"), r.astype("d"))
    gr.SetLineWidth(3); gr.SetLineColor(ROOT.kBlack); gr.Draw("L SAME")
    line = ROOT.TLine(0, 1, min(xmax, 40), 1); line.SetLineStyle(2); line.Draw()
    c.Print(pdf)
c.Print(pdf + "]")

np.set_printoptions(linewidth=200)
for s in summary:
    print(f"\n=== {s['name']} ===  clusters={s['total']}  fail(ps<0|nbins=0)={s['fail']}  E<{EMIN}={s['lowE']}")
    print(f"  kept(ps>=0)={s['n_ge']}  ps>0={s['n_gt']}  ps==0={s['n_eq']}  -> ps==0 fraction {s['n_eq']/max(s['n_ge'],1):.1%}")
    print("  ps histogram (0,0.5,1.0,...):", s["ps_hist"])
    print(f"  {'sel':7s} {'peak[X0]':>8s} {'mean[X0]':>8s} {'<E>':>7s} {'<disc>':>7s} {'<innerOff>':>10s} {'bin0':>7s} {'fγ':>5s}")
    for lab in ("ps>0", "ps>=0", "ps==0"):
        if lab in s:
            d = s[lab]
            print(f"  {lab:7s} {d['peak']:8.2f} {d['mean']:8.2f} {d['E']:7.2f} {d['disc']:7.3f} {d['inner']:10.2f} {d['first_bin']:7.4f} {d['photon']:5.2f}")
print(f"\nPlots: {pdf}")
