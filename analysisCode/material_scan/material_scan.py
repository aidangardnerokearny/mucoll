#!/usr/bin/env python3
"""Material budget vs polar angle from a geantino scan -> a ROOT file of histograms.

Reads the output of a geantino run made with g4sim's material-budget branches
(step_nX0, step_nLambda) and writes one ROOT file. Each step is assigned to a
subdetector after the fact by looking up its midpoint in the GDML geometry
(TGeo), so the run does not need to store volume names:

  rays                 TH1D      geantinos per theta bin
  total_X0             TProfile  X/X0 per ray vs theta (error = spread in bin)
  total_lambda         TProfile  lambda/lambda_I per ray vs theta
  X0/<subdet>          TH1D      mean X/X0 per ray vs theta, one per subdetector
  lambda/<subdet>      TH1D      mean lambda/lambda_I per ray vs theta

Subdetectors are the top-level volumes (daughters of the world) with
"_envelope" and any placement number stripped (ECalBarrel_envelope_11 -> ECalBarrel); the nozzle pieces
are summed into "Nozzle". The per-subdetector histograms add up to the totals,
so e.g. X0/ECalBarrel + X0/ECalEndcap is the whole ECal.

Each event is one straight ray from the gun position, so theta comes from the
primary's initial momentum. Sums run from the origin to where the geantino
leaves the world; --rmax/--zmax stop earlier (e.g. at the calorimeter face).
The step loop runs in C++ (TTreeReader, JIT-compiled by cling).

    python material_scan.py geantino_scan/output_0.root --gdml geantino_scan/MAIA_v0.gdml
    python material_scan.py out.root --gdml MAIA_v0.gdml -o mat.root --rmax 1500 --zmax 2300
    python material_scan.py out.root --gdml MAIA_v0.gdml --exclude Nozzle
"""
import argparse
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
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"

namespace matscan {

// 'ECalBarrel_envelope_11' -> 'ECalBarrel', 'NozzleW_left_2' -> 'Nozzle',
// 'OuterTrackers_9' -> 'OuterTrackers'.
std::string Subdet(std::string v) {
    auto us = v.rfind('_');
    if (us != std::string::npos && us + 1 < v.size()) {
        bool digits = true;
        for (size_t i = us + 1; i < v.size(); ++i)
            if (!std::isdigit(static_cast<unsigned char>(v[i]))) { digits = false; break; }
        if (digits) v.erase(us);
    }
    auto p = v.find("_envelope");
    if (p != std::string::npos) v.erase(p, 9);
    if (v.rfind("Nozzle", 0) == 0) return "Nozzle";
    return v;
}

// Name of the top-level volume (daughter of the world) containing a point in mm,
// e.g. "ECalBarrel_envelope". TGeo works in cm. Level 0 is the world itself.
std::string TopVolume(TGeoManager* geo, double x, double y, double z) {
    geo->FindNode(x / 10., y / 10., z / 10.);
    const int level = geo->GetLevel();
    if (level == 0) return geo->GetTopVolume()->GetName();
    return geo->GetMother(level - 1)->GetVolume()->GetName();
}

struct Result {
    TH1D* rays = nullptr;
    TProfile* totX0 = nullptr;
    TProfile* totLam = nullptr;
    std::map<std::string, TH1D*> x0, lam;
    Long64_t nEvents = 0;
};

const char* kX0Axis  = ";#theta [deg];mean X/X_{0} per ray";
const char* kLamAxis = ";#theta [deg];mean #lambda/#lambda_{I} per ray";

TH1D* Book(const std::string& name, const std::string& title, int nbins) {
    auto* h = new TH1D(name.c_str(), title.c_str(), nbins, 0., 180.);
    h->SetDirectory(nullptr);
    return h;
}

Result Fill(TTree* tree, TGeoManager* geo, int nbins, double rmax, double zmax,
            Long64_t maxEvents, const std::set<std::string>& exclude) {
    Result r;
    r.rays = Book("rays", "geantinos per bin;#theta [deg];rays", nbins);
    r.totX0 = new TProfile("total_X0", (std::string("total X/X_{0}") + kX0Axis).c_str(),
                           nbins, 0., 180., "s");
    r.totLam = new TProfile("total_lambda", (std::string("total #lambda/#lambda_{I}") + kLamAxis).c_str(),
                            nbins, 0., 180., "s");
    r.totX0->SetDirectory(nullptr);
    r.totLam->SetDirectory(nullptr);

    TTreeReader rd(tree);
    TTreeReaderValue<double> px(rd, "primaryStartPx"), py(rd, "primaryStartPy"),
                             pz(rd, "primaryStartPz");
    TTreeReaderValue<std::vector<int>> tid(rd, "step_trackID");
    TTreeReaderValue<std::vector<double>> sx(rd, "step_x"), sy(rd, "step_y"),
                                          sz(rd, "step_z"), slen(rd, "step_length"),
                                          nx0(rd, "step_nX0"), nla(rd, "step_nLambda");

    std::map<std::string, std::string> subdetOf;        // cache volume -> subdet
    std::map<std::string, std::pair<double, double>> perRay;
    while (rd.Next()) {
        if (maxEvents >= 0 && r.nEvents >= maxEvents) break;
        ++r.nEvents;
        const double p = std::sqrt((*px) * (*px) + (*py) * (*py) + (*pz) * (*pz));
        if (p <= 0) continue;
        const double theta = std::acos(std::max(-1., std::min(1., *pz / p))) * 180. / M_PI;
        r.rays->Fill(theta);
        const double ux = *px / p, uy = *py / p, uz = *pz / p;    // straight-line ray

        perRay.clear();
        double sumX0 = 0, sumLam = 0;
        const size_t n = tid->size();
        for (size_t i = 0; i < n; ++i) {
            if ((*tid)[i] != 1) continue;                    // the geantino itself
            if (rmax > 0 && std::hypot((*sx)[i], (*sy)[i]) >= rmax) continue;
            if (zmax > 0 && std::fabs((*sz)[i]) >= zmax) continue;
            // The stored point is the pre-step point, which for a geantino sits ON a
            // boundary; look up the step's midpoint so it lands inside the volume.
            const double h = 0.5 * (*slen)[i];
            const std::string v = TopVolume(geo, (*sx)[i] + h * ux, (*sy)[i] + h * uy,
                                            (*sz)[i] + h * uz);
            auto it = subdetOf.find(v);
            if (it == subdetOf.end()) it = subdetOf.emplace(v, Subdet(v)).first;
            if (exclude.count(it->second)) continue;
            auto& acc = perRay[it->second];
            acc.first += (*nx0)[i];
            acc.second += (*nla)[i];
            sumX0 += (*nx0)[i];
            sumLam += (*nla)[i];
        }
        for (auto& [det, acc] : perRay) {
            auto hx = r.x0.find(det);
            if (hx == r.x0.end()) {
                hx = r.x0.emplace(det, Book(det, det + kX0Axis, nbins)).first;
                r.lam.emplace(det, Book(det, det + kLamAxis, nbins));
            }
            hx->second->Fill(theta, acc.first);
            r.lam[det]->Fill(theta, acc.second);
        }
        r.totX0->Fill(theta, sumX0);
        r.totLam->Fill(theta, sumLam);
    }
    // sums -> mean per ray in each theta bin (rays that miss a subdetector count as 0)
    for (auto& [det, h] : r.x0) { h->Divide(r.rays); r.lam[det]->Divide(r.rays); }
    return r;
}

}  // namespace matscan
""")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", help="g4sim output (output_0.root)")
    ap.add_argument("--gdml", required=True,
                    help="the SAME geometry the run used (e.g. geantino_scan/MAIA_v0.gdml)")
    ap.add_argument("-o", "--out", help="output ROOT file (default: <input>_material.root)")
    ap.add_argument("--tree", default="tree")
    ap.add_argument("--bins", type=int, default=180, help="theta bins over 0-180 deg")
    ap.add_argument("--rmax", type=float, help="only count steps with rho < RMAX [mm]")
    ap.add_argument("--zmax", type=float, help="only count steps with |z| < ZMAX [mm]")
    ap.add_argument("--max-events", type=int, default=-1)
    ap.add_argument("--exclude", nargs="+", default=[], metavar="SUBDET",
                    help="subdetectors to leave out of everything, totals included "
                         "(e.g. --exclude Nozzle)")
    args = ap.parse_args()

    fin = ROOT.TFile.Open(args.root)
    if not fin or fin.IsZombie():
        raise SystemExit(f"cannot open {args.root}")
    tree = fin.Get(args.tree)
    if not tree:
        raise SystemExit(f"{args.root}: no TTree '{args.tree}'")
    need = ["primaryStartPx", "primaryStartPy", "primaryStartPz", "step_trackID",
            "step_x", "step_y", "step_z", "step_length", "step_nX0", "step_nLambda"]
    missing = [b for b in need if not tree.GetBranch(b)]
    if missing:
        raise SystemExit(f"{args.root}: missing branch(es) {', '.join(missing)} -- was it "
                         "written by a g4sim with the material-budget patch?")

    ROOT.TGeoManager.SetVerboseLevel(0)
    geo = ROOT.TGeoManager.Import(args.gdml)
    if not geo:
        raise SystemExit(f"could not load geometry from {args.gdml}")

    exclude = ROOT.std.set("std::string")()
    for d in args.exclude:
        exclude.insert(d)
    res = ROOT.matscan.Fill(tree, geo, args.bins, args.rmax or -1.0, args.zmax or -1.0,
                            args.max_events, exclude)

    out = args.out or str(Path(args.root).with_name(Path(args.root).stem + "_material.root"))
    fout = ROOT.TFile(out, "RECREATE")
    for h in (res.rays, res.totX0, res.totLam):
        h.Write()
    for dirname, hists in (("X0", res.x0), ("lambda", res.lam)):
        d = fout.mkdir(dirname)
        d.cd()
        for name, h in hists:
            h.Write()
        fout.cd()
    fout.Close()

    dets = sorted((str(k) for k, _ in res.x0), key=lambda d: -res.x0[d].Integral())
    print(f"{res.nEvents} rays -> {out}")
    print(f"  {'subdetector':20s} {'max X/X0':>9s} {'max l/lI':>9s}")
    for d in dets:
        print(f"  {d:20s} {res.x0[d].GetMaximum():9.2f} {res.lam[d].GetMaximum():9.2f}")


if __name__ == "__main__":
    main()
