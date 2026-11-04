#!/usr/bin/env python3
"""Plot the histogram file written by material_scan.py.

Two stacked panels vs theta: mean X/X0 per ray (top) and lambda/lambda_I (bottom),
one band per subdetector, with the per-ray total (TProfile) drawn on top.
Barrel and endcap histograms are added together (ECalBarrel + ECalEndcap -> ECal).

    python plot_material.py output_0_material.root                 # -> output_0_material.pdf/.png
    python plot_material.py output_0_material.root --exclude Nozzle --ymax 250 20
"""
import argparse
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# fixed colour order, one per band (never cycled)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
          "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def group(name):
    """ECalBarrel / ECalEndcap -> ECal; everything else unchanged."""
    return name.replace("Barrel", "").replace("Endcap", "")


def load(f, dirname, exclude):
    """{group: TH1D} from one directory of the file, barrel+endcap summed."""
    out = {}
    for key in f.Get(dirname).GetListOfKeys():
        h = key.ReadObj()
        g = group(h.GetName())
        if g in exclude or h.Integral() == 0:
            continue
        if g in out:
            out[g].Add(h)
        else:
            out[g] = h.Clone(f"{dirname}_{g}")
            out[g].SetDirectory(0)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="file from material_scan.py")
    ap.add_argument("-o", "--out", help="output name without extension")
    ap.add_argument("--exclude", nargs="+", default=[], help="groups to leave out, e.g. Nozzle")
    ap.add_argument("--ymax", type=float, nargs=2, metavar=("X0", "LAMBDA"))
    args = ap.parse_args()

    f = ROOT.TFile.Open(args.root)
    out = args.out or str(Path(args.root).with_suffix(""))

    c = ROOT.TCanvas("c", "", 1000, 900)
    c.Divide(1, 2, 0, 0)
    keep = []                                   # PyROOT: keep drawn objects alive
    legend = ROOT.TLegend(0.12, 0.74, 0.90, 0.93)
    legend.SetNColumns(4)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    panels = [("X0", "mean X/X_{0} per ray", "total_X0"),
              ("lambda", "mean #lambda/#lambda_{I} per ray", "total_lambda")]

    # one colour per subdetector, fixed by its X0 rank, used in BOTH panels
    x0 = load(f, "X0", args.exclude)
    lam = load(f, "lambda", args.exclude)
    order = sorted(x0, key=lambda g: -x0[g].Integral())          # biggest at the bottom
    if len(order) > len(COLORS):
        keep, rest = order[:len(COLORS) - 1], order[len(COLORS) - 1:]
        for hists in (x0, lam):                                    # same folding in both panels
            other = hists[rest[0]].Clone(f"Other_{id(hists)}")
            other.SetDirectory(0)
            for g in rest[1:]:
                if g in hists:
                    other.Add(hists[g])
            for g in rest:
                hists.pop(g, None)
            hists["Other"] = other
        order = keep + ["Other"]
        print("folded into Other:", ", ".join(rest))
    color_of = dict(zip(order, COLORS))

    for i, (dirname, ytitle, totname) in enumerate(panels, start=1):
        pad = c.cd(i)
        pad.SetLeftMargin(0.10)
        pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.06 if i == 1 else 0.02)
        pad.SetBottomMargin(0.02 if i == 1 else 0.14)
        #pad.SetGridy()

        hists = x0 if dirname == "X0" else load(f, dirname, args.exclude)

        stack = ROOT.THStack(f"stack_{dirname}", "")
        for g in order:
            if g not in hists:
                continue
            h = hists[g]
            h.SetFillColor(ROOT.TColor.GetColor(color_of[g]))
            h.SetLineColor(ROOT.kWhite)
            stack.Add(h, "hist")
            if i == 1:
                legend.AddEntry(h, g, "f")
        stack.Draw()

        # total with its spread; only meaningful as a cross-check without --exclude
        total = f.Get(totname)
        if not args.exclude:
            total.SetMarkerStyle(20)
            total.SetMarkerSize(0.4)
            total.SetMarkerColor(ROOT.kBlack)
            total.SetLineColor(ROOT.kBlack)
            total.Draw("same e1")
            if i == 1:
                legend.AddEntry(total, "total (spread in bin)", "pe")

        top = args.ymax[i - 1] if args.ymax else stack.GetMaximum() * (1.45 if i == 1 else 1.1)
        stack.SetMaximum(top)
        stack.GetYaxis().SetTitle(ytitle)
        stack.GetYaxis().SetTitleSize(0.055)
        stack.GetYaxis().SetTitleOffset(0.8)
        stack.GetYaxis().SetLabelSize(0.045)
        stack.GetXaxis().SetTitle("#theta [deg]" if i == 2 else "")
        stack.GetXaxis().SetTitleSize(0.06)
        stack.GetXaxis().SetLabelSize(0.05 if i == 2 else 0)
        if i == 1:
            legend.Draw()
        pad.RedrawAxis()
        keep += [stack, hists, total]

    for ext in ("pdf", "png"):
        c.SaveAs(f"{out}.{ext}")


if __name__ == "__main__":
    main()
