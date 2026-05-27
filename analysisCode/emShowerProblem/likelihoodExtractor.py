"""
plot_pandora_likelihoods.py
----------------------------
Parse and plot the photon ID likelihood PDFs stored in
PandoraLikelihoodData12EBin.xml (from madbaron/SteeringMacros).

Usage
-----
  python plot_pandora_likelihoods.py [path/to/PandoraLikelihoodData12EBin.xml]

If no path is given it looks for the file in the current directory.
Plots are saved as PNGs alongside the script, and also shown interactively
if a display is available.

Output files
------------
  pandora_<variable>_all_ebins.png   — one canvas per variable, all energy bins
  pandora_all_vars_ebin<N>.png       — one canvas per energy bin, all variables
  pandora_likelihood_ratios.png      — log(S/B) ratio for every variable × energy bin
"""

import sys
import os
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless-safe; remove this line if you want interactive windows
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# ── pretty names and axis labels ─────────────────────────────────────────────

VARIABLE_LABELS = {
    "PeakRms":                ("Peak RMS",               r"Transverse RMS at shower peak  [cm]"),
    "LongProfileStart":       ("Long. Profile Start",    r"Shower start layer"),
    "LongProfileDiscrepancy": ("Long. Profile Discr.",   r"Longitudinal profile discrepancy"),
    "PeakEnergyFraction":     ("Peak Energy Fraction",   r"Fraction of energy in shower peak"),
    "RmsRatio":               ("RMS Ratio",              r"Peak RMS / overall cluster RMS"),
    "MinDistanceToTrack":     ("Min Track Distance",     r"Minimum distance to nearest track  [mm]"),
}

# ── parsing ───────────────────────────────────────────────────────────────────

def load_xml(path):
    with open(path) as f:
        raw = f.read()
    # File has no single root element — wrap it
    return ET.fromstring(f"<root>{raw}</root>")


def parse_global(root):
    n_ebins     = int(root.find("NEnergyBins").text)
    e_edges     = list(map(float, root.find("EnergyBinLowerEdges").text.split()))
    n_sig       = list(map(int,   root.find("NSignalEvents").text.split()))
    n_bkg       = list(map(int,   root.find("NBackgroundEvents").text.split()))
    return n_ebins, e_edges, n_sig, n_bkg


def parse_histogram(root, tag):
    """Return (bin_edges, bin_contents) for the named tag, or None."""
    el = root.find(tag)
    if el is None:
        return None
    n   = int(el.find("NBinsX").text)
    lo  = float(el.find("XLow").text)
    hi  = float(el.find("XHigh").text)
    raw = list(map(float, el.find("BinContents").text.split()))
    # The XML sometimes stores n+1 values (overflow sentinel); keep only n
    contents = np.array(raw[:n])
    edges    = np.linspace(lo, hi, n + 1)
    return edges, contents


def load_all_histograms(root, variables, n_ebins):
    """
    Returns a nested dict:
      data[variable][ebin] = {"sig": (edges, contents), "bkg": (edges, contents)}
    """
    data = {}
    for var in variables:
        data[var] = {}
        for i in range(n_ebins):
            sig = parse_histogram(root, f"PhotonSig{var}_{i}")
            bkg = parse_histogram(root, f"PhotonBkg{var}_{i}")
            if sig and bkg:
                data[var][i] = {"sig": sig, "bkg": bkg}
    return data


# ── helpers ───────────────────────────────────────────────────────────────────

def ebin_label(i, e_edges, n_ebins):
    lo = e_edges[i]
    hi = e_edges[i + 1] if i + 1 < len(e_edges) else "∞"
    if isinstance(hi, float):
        return f"{lo}–{hi} GeV"
    return f">{lo} GeV"


def draw_hist(ax, edges, contents, label, color, alpha=0.75, lw=1.5):
    centres = 0.5 * (edges[:-1] + edges[1:])
    ax.step(edges[:-1], contents, where="post", color=color,
            linewidth=lw, label=label, alpha=alpha)
    ax.fill_between(edges[:-1], contents, step="post",
                    color=color, alpha=alpha * 0.25)


# ── plot 1: one figure per variable, all energy bins ─────────────────────────

def plot_per_variable(data, variables, n_ebins, e_edges, out_dir):
    cmap = cm.get_cmap("plasma", n_ebins)

    for var in variables:
        nice, xlabel = VARIABLE_LABELS.get(var, (var, var))
        ncols = 4
        nrows = (n_ebins + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(ncols * 3.8, nrows * 3.0),
                                 constrained_layout=True)
        fig.suptitle(f"Photon ID PDFs — {nice}", fontsize=14, fontweight="bold")
        axes_flat = axes.flatten()

        for i in range(n_ebins):
            ax = axes_flat[i]
            color = cmap(i / max(n_ebins - 1, 1))
            hists = data[var].get(i)
            if hists:
                draw_hist(ax, *hists["sig"], label="Signal (γ)", color="#2196F3")
                draw_hist(ax, *hists["bkg"], label="Background",  color="#F44336")
            ax.set_title(ebin_label(i, e_edges, n_ebins), fontsize=8)
            ax.set_xlabel(xlabel, fontsize=7)
            ax.set_ylabel("Normalised density", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.set_ylim(bottom=0)
            if i == 0:
                ax.legend(fontsize=6, framealpha=0.5)

        # hide unused panels
        for j in range(n_ebins, len(axes_flat)):
            axes_flat[j].set_visible(False)

        fname = os.path.join(out_dir, f"pandora_{var}_all_ebins.png")
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"  Saved {fname}")


# ── plot 2: one figure per energy bin, all variables ─────────────────────────

def plot_per_ebin(data, variables, n_ebins, e_edges, out_dir):
    for i in range(n_ebins):
        nv = len(variables)
        ncols = 3
        nrows = (nv + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(ncols * 4.0, nrows * 3.2),
                                 constrained_layout=True)
        fig.suptitle(f"Photon ID PDFs — {ebin_label(i, e_edges, n_ebins)}",
                     fontsize=13, fontweight="bold")
        axes_flat = axes.flatten()

        for j, var in enumerate(variables):
            ax = axes_flat[j]
            nice, xlabel = VARIABLE_LABELS.get(var, (var, var))
            hists = data[var].get(i)
            if hists:
                draw_hist(ax, *hists["sig"], label="Signal (γ)", color="#2196F3")
                draw_hist(ax, *hists["bkg"], label="Background",  color="#F44336")
            ax.set_title(nice, fontsize=9, fontweight="bold")
            ax.set_xlabel(xlabel, fontsize=8)
            ax.set_ylabel("Normalised density", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.set_ylim(bottom=0)
            ax.legend(fontsize=7, framealpha=0.5)

        for j in range(nv, len(axes_flat)):
            axes_flat[j].set_visible(False)

        fname = os.path.join(out_dir, f"pandora_all_vars_ebin{i:02d}.png")
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"  Saved {fname}")


# ── plot 3: log(S/B) likelihood ratio heat overview ──────────────────────────

def plot_likelihood_ratios(data, variables, n_ebins, e_edges, out_dir):
    """
    For each (variable, energy bin) pair, overlay log(S/B) as a function of
    the observable value.  Displayed as a grid of line plots, coloured by
    energy bin.
    """
    nv   = len(variables)
    ncols = 3
    nrows = (nv + ncols - 1) // ncols
    cmap  = cm.get_cmap("viridis", n_ebins)

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 4.5, nrows * 3.5),
                             constrained_layout=True)
    fig.suptitle("Log-likelihood ratio  log(S/B)  per variable & energy bin",
                 fontsize=13, fontweight="bold")
    axes_flat = axes.flatten()

    sm = cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=n_ebins - 1))
    sm.set_array([])

    for j, var in enumerate(variables):
        ax = axes_flat[j]
        nice, xlabel = VARIABLE_LABELS.get(var, (var, var))

        for i in range(n_ebins):
            hists = data[var].get(i)
            if not hists:
                continue
            edges, sig = hists["sig"]
            _,     bkg = hists["bkg"]
            centres = 0.5 * (edges[:-1] + edges[1:])
            eps = 1e-9
            ratio = np.log((sig + eps) / (bkg + eps))
            ax.plot(centres, ratio, color=cmap(i / max(n_ebins - 1, 1)),
                    linewidth=1.2, alpha=0.85,
                    label=ebin_label(i, e_edges, n_ebins))

        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax.set_title(nice, fontsize=9, fontweight="bold")
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("log(S / B)", fontsize=8)
        ax.tick_params(labelsize=7)

    for j in range(nv, len(axes_flat)):
        axes_flat[j].set_visible(False)

    cbar = fig.colorbar(sm, ax=axes_flat[:nv], shrink=0.6, pad=0.02,
                        label="Energy bin index")
    # annotate a few energy labels on the colorbar
    cbar.set_ticks(range(n_ebins))
    cbar.set_ticklabels([ebin_label(i, e_edges, n_ebins) for i in range(n_ebins)],
                        fontsize=5)

    fname = os.path.join(out_dir, "pandora_likelihood_ratios.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  Saved {fname}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    xml_path = sys.argv[1] if len(sys.argv) > 1 else "PandoraLikelihoodData12EBin.xml"

    if not os.path.exists(xml_path):
        print(f"ERROR: Cannot find {xml_path}")
        print("Usage: python plot_pandora_likelihoods.py [path/to/PandoraLikelihoodData12EBin.xml]")
        sys.exit(1)

    #out_dir = os.path.dirname(os.path.abspath(xml_path))
    out_dir = "/home/aidangardnerokearny/mucoll/analysisCode/emShowerProblem/showerProblemPlots/likelihoodPlots"
    print(f"Reading {xml_path} ...")

    root              = load_xml(xml_path)
    n_ebins, e_edges, n_sig, n_bkg = parse_global(root)

    print(f"  {n_ebins} energy bins: {e_edges}")
    print(f"  Signal events per bin:     {n_sig}")
    print(f"  Background events per bin: {n_bkg}")

    # Discover variables from tag names present in the file
    variables_found = set()
    for child in root:
        tag = child.tag
        for prefix in ("PhotonSig", "PhotonBkg"):
            if tag.startswith(prefix):
                core = tag[len(prefix):]              # e.g. "PeakRms_0"
                varname = "_".join(core.split("_")[:-1])   # strip trailing _N
                if varname:
                    variables_found.add(varname)

    # Canonical ordering
    canonical_order = ["PeakRms", "LongProfileStart", "LongProfileDiscrepancy",
                       "PeakEnergyFraction", "RmsRatio", "MinDistanceToTrack"]
    variables = [v for v in canonical_order if v in variables_found]
    variables += sorted(variables_found - set(canonical_order))   # any extras
    print(f"  Variables found: {variables}")

    data = load_all_histograms(root, variables, n_ebins)

    print("\nGenerating plots ...")
    plot_per_variable(data, variables, n_ebins, e_edges, out_dir)
    plot_per_ebin(data, variables, n_ebins, e_edges, out_dir)
    plot_likelihood_ratios(data, variables, n_ebins, e_edges, out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
