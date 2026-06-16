"""
pidVarPlotter.jl

Plots the 4 shower-profile classifier variables available in the tree,
split by whether the PFO was tagged as a photon (pdgCode == 22) or
something else. Mirrors plot_pid_separation.py exactly.

Variables plotted:
  1. profileStart       (LongProfileStart in likelihood XML)
  2. profileDiscrepancy (LongProfileDiscrepancy in likelihood XML)
  3. clusterEnergy      (bonus: energy scale context)
  4. showerStartOffset  (showerStartLayer - innerLayer: pre-showering proxy)
  5. 2D contour: profileStart vs profileDiscrepancy

All 1D distributions are normalised to unit area.

Usage
-----
    julia plot_pid_separation.jl [input.root] [output.pdf] [e_min] [e_max]

    Defaults: shower_profiles.root  pid_separation.pdf  0.5  1e9

Dependencies (install once):
    using Pkg
    Pkg.add(["UnROOT", "Plots", "StatsBase"])
"""

using UnROOT
using Plots
using StatsBase
using Statistics
using Printf

# ── Arguments ─────────────────────────────────────────────────────────────────
input_file  = length(ARGS) >= 1 ? ARGS[1] : "shower_profiles.root"
output_file = length(ARGS) >= 2 ? ARGS[2] : "pid_separation.pdf"
energy_min  = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.5
energy_max  = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1e9

# ── Global plot defaults ───────────────────────────────────────────────────────
default(
    grid       = true,
    gridalpha  = 0.3,
    framestyle = :box,
    legend_background_color = :transparent,
    legend_foreground_color = :transparent,
    size       = (900, 600),
)

# ── Open ROOT file ─────────────────────────────────────────────────────────────
println("Opening $input_file ...")
f    = ROOTFile(input_file)

available = keys(f["ShowerProfiles"])
println(" Available branches: ", join(available, ", "))

wanted =  [
    "clusterEnergy", "profileStart", "profileDiscrepancy",
    "showerStartOffset", "innerLayer", "showerStartLayer",
    "nProfileBins", "pdgCode", "peakRms", "rmsRatio"
]
present = filter(b->b in available, wanted)
missing_branches = filter(b-> !(b in available), wanted)

isempty(missing_branches) || println(" WARNING: branches not found and will be skipped: ", join(missing_branches, ", "))

tree = LazyTree(f, "ShowerProfiles", present)
println("  $(length(tree)) entries")

has_peakrms  = "peakRms"          in present
has_rrms = "rmsRatio"         in present
has_pdgcode  = "pdgCode"          in present
has_sso      = "showerStartLayer" in present

# ── Collect data split by PID ──────────────────────────────────────────────────
ph_ps   = Float64[]; oth_ps   = Float64[]
ph_pd   = Float64[]; oth_pd   = Float64[]
ph_e    = Float64[]; oth_e    = Float64[]
ph_sso  = Float64[]; oth_sso  = Float64[]
ph_prms = Float64[]; oth_prms = Float64[]
ph_rrms = Float64[]; oth_rrms = Float64[]

n_no_pdg = 0

for row in tree
    E = Float64(row.clusterEnergy)
    (E < energy_min || E > energy_max) && continue
    row.profileStart <= 0               && continue
    row.nProfileBins == 0               && continue

    # pdgCode is 0 if branch was missing (UnROOT returns default)
    pdg = Int(row.pdgCode)
    if pdg == 0
        n_no_pdg += 1
    end

    is_photon = pdg == 22

    ps   = Float64(row.profileStart)
    pd   = Float64(row.profileDiscrepancy)
    sso  = Float64(row.showerStartOffset)
    prms = has_peakrms ? Float64(row.peakRms) : -1.0
    rrms = has_rrms ? Float64(row.rmsRatio)   : -1.0

    if is_photon
        push!(ph_ps,  ps);  push!(ph_pd,  pd)
        push!(ph_e,   E);   push!(ph_sso, sso)
        has_peakrms  && push!(ph_prms, prms)
        has_rrms && push!(ph_rrms, rrms)
    else
        push!(oth_ps,  ps);  push!(oth_pd,  pd)
        push!(oth_e,   E);   push!(oth_sso, sso)
        has_peakrms  && push!(oth_prms, prms)
        has_rrms && push!(oth_rrms, rrms)
    end
end

n_no_pdg > 0 && println("WARNING: $n_no_pdg entries had pdgCode=0 — assigned to 'other'")
n_ph  = length(ph_ps)
n_oth = length(oth_ps)
println("  Photon-tagged (pdg=22) : $n_ph")
println("  Other                  : $n_oth")
(n_ph == 0 && n_oth == 0) && error("No entries pass selection.")

# ── Helper: normalised overlaid histogram ─────────────────────────────────────
function sep_plot(title, xlabel, ph_vals, oth_vals;
                  nbins=50, x_lo=nothing, x_hi=nothing, logy=false)

    all_vals = [ph_vals; oth_vals]
    isempty(all_vals) && return plot(; title, xlabel)

    xlo = isnothing(x_lo) ? minimum(all_vals) : x_lo
    xhi = isnothing(x_hi) ? maximum(all_vals) * 1.05 : x_hi

    edges = range(xlo, xhi, length = nbins + 1)

    # Normalise to unit area
    function norm_hist(vals)
        h = fit(Histogram, Float64.(vals), edges)
        s = sum(h.weights)
        s > 0 ? (h.weights ./ s) : h.weights
    end

    ph_w  = norm_hist(ph_vals)
    oth_w = norm_hist(oth_vals)

    centres = collect(edges[1:end-1]) .+ step(edges) / 2

    p = plot(;
        title,
        xlabel,
        ylabel = "Normalised entries / bin",
        legend = :topright,
        yscale = logy ? :log10 : :identity,
        ylims  = logy ? (1e-4, :auto) : (0, :auto),
    )

    bar!(p, centres, ph_w;
        bar_width  = step(edges),
        fillalpha  = 0.3,
        fillcolor  = :royalblue,
        linecolor  = :royalblue,
        linewidth  = 1.5,
        label      = "Photon (n=$n_ph)",
    )
    bar!(p, centres, oth_w;
        bar_width  = step(edges),
        fillalpha  = 0.3,
        fillcolor  = :crimson,
        linecolor  = :crimson,
        linewidth  = 1.5,
        label      = "Other (n=$n_oth)",
    )

    return p
end

# ── Helper: 2D contour overlay ────────────────────────────────────────────────
function sep_2d(ph_x, ph_y, oth_x, oth_y; 
        nbins=50, xlabel="x", ylabel="y", title="2D")
    all_x = [ph_x; oth_x]
    all_y = [ph_y; oth_y]

    x_edges = range(0, maximum(all_x) * 1.05, length = nbins + 1)
    y_hi    = min(maximum(all_y) * 1.05, 5.0)
    y_edges = range(0, y_hi, length = nbins + 1)

    function norm_h2(xs, ys)
        h = fit(Histogram, (Float64.(xs), Float64.(ys)),
                (x_edges, y_edges))
        s = sum(h.weights)
        s > 0 ? h.weights ./ s : h.weights
    end

    ph_w  = norm_h2(ph_x,  ph_y)
    oth_w = norm_h2(oth_x, oth_y)

    xc = collect(x_edges[1:end-1]) .+ step(x_edges) / 2
    yc = collect(y_edges[1:end-1]) .+ step(y_edges) / 2

    p = contour(xc, yc, ph_w';
        levels    = 8,
        color     = :blues,
        linewidth = 2,
        label     = "Photon (n=$n_ph)",
        xlabel    = "profileStart (X₀)",
        ylabel    = "profileDiscrepancy (normalised)",
        title     = "profileStart vs profileDiscrepancy",
        colorbar  = false,
    )
    contour!(p, xc, yc, oth_w';
        levels    = 8,
        color     = :reds,
        linewidth = 2,
        label     = "Other (n=$n_oth)",
        colorbar  = false,
    )

    return p
end

# ── Build pages ────────────────────────────────────────────────────────────────
all_ps   = [ph_ps;  oth_ps]
all_pd   = [ph_pd;  oth_pd]
all_e    = [ph_e;   oth_e]
all_sso  = [ph_sso; oth_sso]
all_prms = [ph_prms; oth_prms]
all_rrms = [ph_rrms; oth_rrms]

pages = Plots.Plot[]

push!(pages, sep_plot(
    "LongProfileStart: photon vs other",
    "profileStart (radiation lengths X₀)",
    ph_ps, oth_ps;
    nbins = 50, x_lo = 0.0, x_hi = maximum(all_ps) * 1.05,
))

push!(pages, sep_plot(
    "LongProfileDiscrepancy: photon vs other",
    "profileDiscrepancy (normalised)",
    ph_pd, oth_pd;
    nbins = 50, x_lo = 0.0, x_hi = min(maximum(all_pd) * 1.05, 5.0),
))

if has_peakrms
    push!(pages, sep_plot(
        "PeakRms: photon vs other",
        "peakRms",
        ph_prms, oth_prms;
        nbins = 50, x_lo = 0.0,
        x_hi = isempty(all_prms) ? 5.0 : min(maximum(filter(v->v>=0,all_prms)) * 1.05, 10.0),
    ))
end

if has_rrms
    push!(pages, sep_plot(
        "RmsRatio (RmsXYRatio): photon vs other",
        "rmsRatio",
        ph_rrms, oth_rrms;
        nbins = 50, x_lo = 0.0,
        x_hi = isempty(all_rrms) ? 5.0 : min(maximum(filter(v->v>=0,all_rrms)) * 1.05, 10.0),
    ))
end

push!(pages, sep_plot(
    "Cluster EM energy: photon vs other",
    "Cluster EM energy [GeV]",
    ph_e, oth_e;
    nbins = 60, x_lo = 0.0, x_hi = maximum(all_e) * 1.05, logy = true,
))

if has_sso
    push!(pages, sep_plot(
        "Shower start offset: photon vs other  (showerStartLayer − innerLayer)",
        "Pseudo-layer offset",
        ph_sso, oth_sso;
        nbins = 60,
        x_lo = Float64(minimum(all_sso)) - 1,
        x_hi = Float64(maximum(all_sso)) + 1,
    ))
end

push!(pages, sep_2d(
    ph_ps, ph_pd, oth_ps, oth_pd;
    xlabel = "profileStart (X₀)",
    ylabel = "profileDiscrepancy (normalised)",
    title  = "profileStart vs profileDiscrepancy",
))

if has_peakrms && has_rrms
    push!(pages, sep_2d(
        ph_prms, ph_rrms, oth_prms, oth_rrms;
        xlabel = "peakRms",
        ylabel = "rmsRatio",
        title  = "peakRms vs rmsRatio",
    ))
end

println(pages)

# ── Save to PDF ────────────────────────────────────────────────────────────────
println("Saving $(length(pages)) pages to $output_file ...")
savefig(pages[1], output_file)
for p in pages[2:end]
    savefig(p, output_file)
end
println("Saved $output_file")
