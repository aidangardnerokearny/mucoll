#!/usr/bin/env python3
"""Material budget vs polar angle from a geantino scan (PyROOT version).

Reads the output of a geantino run (e.g. maia_geantino_theta_scan.yaml) made
with g4sim's material-budget branches (step_nX0, step_nLambda, step_volume)
and writes, next to each other in --outdir:

  material_vs_theta.pdf/.png   mean X/X0 (top) and lambda/lambda_I (bottom) per
                               ray vs theta, stacked by subdetector
  material_vs_theta.root       the histograms behind the plot (TH1D per
                               subdetector + a TProfile of the per-ray total)
  material_vs_theta.csv        the same numbers as a table

Each event is one straight ray from the gun position, so theta comes from the
primary's initial momentum. Sums run from the origin to where the geantino
leaves the world; --rmax/--zmax stop earlier (e.g. at the calorimeter face).
The step loop runs in C++ (TTreeReader, JIT-compiled by cling), so large scans
are fast.

    python material_scan.py geantino_scan/output_0.root -o geantino_scan/material
    python material_scan.py out.root --bins 90 --rmax 1500 --zmax 2300   # mm
    python material_scan.py out.root --exclude Nozzle      # zoom on the detector
"""
import argparse
import csv
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="\nThis distribution of ROOT")  # pip-wheel ROOT
import ROOT  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

ROOT.gInterpreter.Declare(r"""
#include <map>
#include <set>
#include <string>
#include <vector>
#include <cmath>
#include <cctype>
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1D.h"
#include "TProfile.h"

namespace matscan {

// 'ECalBarrel_envelope_11' -> 'ECal', 'NozzleW_left_2' -> 'Nozzle',
// 'OuterTrackers_9' -> 'OuterTrackers'. Barrel and endcap merge: theta
// already separates them, and it keeps the stack to <= 8 colours.
std::string Subdet(std::string v) {
    auto us = v.rfind('_');
    if (us != std::string::npos && us + 1 < v.size()) {
        bool digits = true;
        for (size_t i = us + 1; i < v.size(); ++i)
            if (!std::isdigit(static_cast<unsigned char>(v[i]))) { digits = false; break; }
        if (digits) v.erase(us);
    }
    for (const char* s : {"_envelope", "Barrel", "Endcap"}) {
        auto p = v.find(s);
        if (p != std::string::npos) v.erase(p, std::char_traits<char>::length(s));
    }
    if (v.rfind("Nozzle", 0) == 0) return "Nozzle";
    return v;
}

struct Result {
    TH1D* rays = nullptr;                       // rays per theta bin
    TProfile* totX0 = nullptr;                  // per-ray total, mean +- spread
    TProfile* totLam = nullptr;
    std::map<std::string, TH1D*> x0, lam;       // summed per subdetector
    Long64_t nEvents = 0;
};

TH1D* Book(const std::string& name, int nbins) {
    auto* h = new TH1D(name.c_str(), "", nbins, 0., 180.);
    h->SetDirectory(nullptr);
    h->Sumw2(false);
    return h;
}

Result Fill(TTree* tree, int nbins, double rmax, double zmax, Long64_t maxEvents,
            const std::set<std::string>& exclude) {
    Result r;
    r.rays = Book("rays", nbins);
    r.totX0 = new TProfile("total_X0", "", nbins, 0., 180., "s");
    r.totLam = new TProfile("total_lambda", "", nbins, 0., 180., "s");
    r.totX0->SetDirectory(nullptr);
    r.totLam->SetDirectory(nullptr);

    TTreeReader rd(tree);
    TTreeReaderValue<double> px(rd, "primaryStartPx"), py(rd, "primaryStartPy"),
                             pz(rd, "primaryStartPz");
    TTreeReaderValue<std::vector<int>> tid(rd, "step_trackID");
    TTreeReaderValue<std::vector<double>> sx(rd, "step_x"), sy(rd, "step_y"),
                                          sz(rd, "step_z"),
                                          nx0(rd, "step_nX0"), nla(rd, "step_nLambda");
    TTreeReaderValue<std::vector<std::string>> vol(rd, "step_volume");

    std::map<std::string, std::string> subdetOf;        // cache volume -> subdet
    std::map<std::string, std::pair<double, double>> perRay;
    while (rd.Next()) {
        if (maxEvents >= 0 && r.nEvents >= maxEvents) break;
        ++r.nEvents;
        const double p = std::sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
        if (p <= 0) continue;
        const double theta = std::acos(std::max(-1., std::min(1., *pz / p))) * 180. / M_PI;
        r.rays->Fill(theta);

        perRay.clear();
        double sumX0 = 0, sumLam = 0;
        const size_t n = tid->size();
        for (size_t i = 0; i < n; ++i) {
            if ((*tid)[i] != 1) continue;                    // the geantino itself
            if (rmax > 0 && std::hypot((*sx)[i], (*sy)[i]) >= rmax) continue;
            if (zmax > 0 && std::fabs((*sz)[i]) >= zmax) continue;
            const std::string& v = (*vol)[i];
            auto it = subdetOf.find(v);
            if (it == subdetOf.end()) it = subdetOf.emplace(v, Subdet(v)).first;
            if (exclude.count(it->second)) continue;         // dropped from everything
            auto& acc = perRay[it->second];
            acc.first += (*nx0)[i];
            acc.second += (*nla)[i];
            sumX0 += (*nx0)[i];
            sumLam += (*nla)[i];
        }
        for (auto& [det, acc] : perRay) {
            auto hx = r.x0.find(det);
            if (hx == r.x0.end()) {
                hx = r.x0.emplace(det, Book("X0_" + det, nbins)).first;
                r.lam.emplace(det, Book("lambda_" + det, nbins));
            }
            hx->second->Fill(theta, acc.first);
            r.lam[det]->Fill(theta, acc.second);
        }
        r.totX0->Fill(theta, sumX0);
        r.totLam->Fill(theta, sumLam);
    }
    // sums -> mean per ray in each theta bin
    for (auto& [det, h] : r.x0) { h->Divide(r.rays); r.lam[det]->Divide(r.rays); }
    return r;
}

}  // namespace matscan
""")

# Fixed categorical order (validated for adjacent-pair CVD separation);
# slot 8 is kept for "Other" only when there are more than 8 subdetectors.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
           "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
OTHER = "#8a8984"
MAX_SERIES = 8


def check_branches(tree, path):
    need = ["primaryStartPx", "primaryStartPy", "primaryStartPz", "step_trackID",
            "step_x", "step_y", "step_z", "step_nX0", "step_nLambda", "step_volume"]
    missing = [b for b in need if not tree.GetBranch(b)]
    if missing:
        raise SystemExit(f"{path}: missing branch(es) {', '.join(missing)} -- was it "
                         "written by a g4sim with the material-budget patch?")


def group_series(res):
    """Order subdetectors by total X0 (largest at the bottom of the stack) and
    fold anything past 7 into 'Other' so hues are never cycled."""
    dets = sorted((str(k) for k, _ in res.x0), key=lambda d: -res.x0[d].Integral())
    dets = [d for d in dets if res.x0[d].Integral() > 0 or res.lam[d].Integral() > 0]
    if len(dets) <= MAX_SERIES:
        return [(d, [d]) for d in dets]
    keep = dets[:MAX_SERIES - 1]
    return [(d, [d]) for d in keep] + [("Other", dets[MAX_SERIES - 1:])]


def summed(hists, members, name):
    h = hists[members[0]].Clone(name)
    h.SetDirectory(0)
    for m in members[1:]:
        h.Add(hists[m])
    return h


def style_stack_member(h, color):
    c = ROOT.TColor.GetColor(color)
    h.SetFillColor(c)
    h.SetLineColor(ROOT.kWhite)       # thin surface gap between stacked segments
    h.SetLineWidth(1)


def draw(series, res, outdir, ymax=None):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetFrameLineColor(ROOT.TColor.GetColor("#6b6a64"))
    ROOT.gStyle.SetGridColor(ROOT.TColor.GetColor("#e4e3dd"))
    ROOT.gStyle.SetGridStyle(1)
    ROOT.gStyle.SetTextFont(42)

    c = ROOT.TCanvas("c", "", 1000, 900)
    c.Divide(1, 2, 0, 0)
    keep = []                                      # keep PyROOT objects alive
    legend = ROOT.TLegend(0.13, 0.72, 0.88, 0.92)
    legend.SetNColumns(4)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.045)

    for pad_i, (key, ytitle, total) in enumerate(
            [("x0", "mean X/X_{0} per ray", res.totX0),
             ("lam", "mean #lambda/#lambda_{I} per ray", res.totLam)], start=1):
        pad = c.cd(pad_i)
        pad.SetLeftMargin(0.1)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.04 if pad_i == 2 else 0.06)
        pad.SetBottomMargin(0.14 if pad_i == 2 else 0.02)
        pad.SetGridy()

        hists = res.x0 if key == "x0" else res.lam
        stack = ROOT.THStack(f"stack_{key}", "")
        for i, (label, members) in enumerate(series):
            h = summed(hists, members, f"{key}_{label}")
            style_stack_member(h, OTHER if label == "Other" else PALETTE[i])
            stack.Add(h, "hist")
            keep.append(h)
            if pad_i == 1:
                legend.AddEntry(h, label, "f")
        stack.Draw()
        if ymax:
            stack.SetMaximum(ymax[pad_i - 1])
        else:   # headroom for the legend, measured against the error bars too
            top = max(stack.GetMaximum(),
                      max(total.GetBinContent(b) + total.GetBinError(b)
                          for b in range(1, total.GetNbinsX() + 1)))
            stack.SetMaximum(top * (1.45 if pad_i == 1 else 1.1))
        ax = stack.GetXaxis()
        ax.SetTitle("#theta [deg]" if pad_i == 2 else "")
        ax.SetLabelSize(0.05 if pad_i == 2 else 0)
        ax.SetTitleSize(0.06)
        ay = stack.GetYaxis()
        ay.SetTitle(ytitle)
        ay.SetTitleSize(0.055)
        ay.SetTitleOffset(0.8)
        ay.SetLabelSize(0.045)

        # per-ray total with its spread inside each theta bin, as a check on the stack
        total.SetMarkerStyle(20)
        total.SetMarkerSize(0.4)
        total.SetMarkerColor(ROOT.TColor.GetColor("#1f1f1e"))
        total.SetLineColor(ROOT.TColor.GetColor("#1f1f1e"))
        total.Draw("same e1")
        if pad_i == 1:
            legend.AddEntry(total, "total (spread in bin)", "pe")
            legend.Draw()
        keep += [stack]
        pad.RedrawAxis()

    for ext in ("pdf", "png"):
        c.SaveAs(str(outdir / f"material_vs_theta.{ext}"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("-o", "--outdir", default="material")
    ap.add_argument("--tree", default="tree")
    ap.add_argument("--bins", type=int, default=180, help="theta bins over 0-180 deg")
    ap.add_argument("--rmax", type=float, help="only count steps with rho < RMAX [mm]")
    ap.add_argument("--zmax", type=float, help="only count steps with |z| < ZMAX [mm]")
    ap.add_argument("--max-events", type=int, default=-1)
    ap.add_argument("--exclude", nargs="+", default=[], metavar="SUBDET",
                    help="subdetectors to leave out of the stack AND the totals, "
                         "e.g. --exclude Nozzle Yoke (names as printed in the summary)")
    ap.add_argument("--ymax", type=float, nargs=2, metavar=("X0", "LAMBDA"),
                    help="fix the y-axis maxima of the two panels")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    f = ROOT.TFile.Open(args.root)
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {args.root}")
    tree = f.Get(args.tree)
    if not tree:
        raise SystemExit(f"{args.root}: no TTree '{args.tree}'")
    check_branches(tree, args.root)

    exclude = ROOT.std.set("std::string")()
    for d in args.exclude:
        exclude.insert(d)
    res = ROOT.matscan.Fill(tree, args.bins, args.rmax or -1.0, args.zmax or -1.0,
                            args.max_events, exclude)
    series = group_series(res)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # histograms, for further work in ROOT
    fout = ROOT.TFile(str(outdir / "material_vs_theta.root"), "RECREATE")
    for h in [res.rays, res.totX0, res.totLam] + \
             [res.x0[str(k)] for k, _ in res.x0] + [res.lam[str(k)] for k, _ in res.lam]:
        h.Write()
    fout.Close()

    # table
    dets = [str(k) for k, _ in res.x0]
    with open(outdir / "material_vs_theta.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["theta_deg", "n_rays", "X0_total", "X0_total_spread",
                    "lambda_total", "lambda_total_spread"]
                   + [f"X0_{d}" for d in dets] + [f"lambda_{d}" for d in dets])
        for b in range(1, args.bins + 1):
            w.writerow([f"{res.rays.GetBinCenter(b):.3f}", int(res.rays.GetBinContent(b)),
                        res.totX0.GetBinContent(b), res.totX0.GetBinError(b),
                        res.totLam.GetBinContent(b), res.totLam.GetBinError(b)]
                       + [res.x0[d].GetBinContent(b) for d in dets]
                       + [res.lam[d].GetBinContent(b) for d in dets])

    b90 = res.totX0.FindBin(90.0)
    print(f"{res.nEvents} rays; X/X0 at theta=90 deg: {res.totX0.GetBinContent(b90):.1f}")
    for d in sorted(dets, key=lambda d: -res.x0[d].Integral()):
        print(f"  {d:20s} max X/X0 {res.x0[d].GetMaximum():8.2f}   "
              f"max l/lI {res.lam[d].GetMaximum():6.2f}")

    if not args.no_plot:
        draw(series, res, outdir, args.ymax)
    print(f"wrote {outdir}/material_vs_theta.{{root,csv}}"
          + ("" if args.no_plot else " and .pdf/.png"))


if __name__ == "__main__":
    main()
