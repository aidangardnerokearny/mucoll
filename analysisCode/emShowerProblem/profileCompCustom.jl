using UnROOT 
using StatsBase
using Statistics
using SpecialFunctions
using LaTeXStrings
using Printf
using CairoMakie

# --- Arguments ---
file_nom = length(ARGS) >= 1 ? ARGS[1] : "shower_profiles.root"
file_vac = length(ARGS) >= 2 ? ARGS[2] : "shower_profiles_vacuum.root"
output_dir = length(ARGS) >= 3 ? ARGS[3] : "/customSVGs"
energy_min = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.1
energy_max = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1e9


# --- Open Files ---
println("Opening $file_nom...")
f_nom = ROOTFile(file_nom)
println("Opening $file_vac...")
f_vac = ROOTFile(file_vac)

tree_nom = LazyTree(file_nom, "ShowerProfiles")
tree_vac = LazyTree(file_vac, "ShowerProfiles")
trees = [tree_nom, tree_vac]

println("Nominal File has $(length(tree_nom)) entries")
println("Vacuum File has $(length(tree_vac)) entries")


# --- Set Up Data/Read From Trees ---

nom_ps     = Float64[]; vac_ps     = Float64[]
nom_pd     = Float64[]; vac_pd     = Float64[]
nom_e      = Float64[]; vac_e      = Float64[]
nom_sso    = Float64[]; vac_sso    = Float64[]
nom_prms   = Float64[]; vac_prms   = Float64[]
nom_rrms   = Float64[]; vac_rrms   = Float64[]
nom_m_prof = Float64[]; vac_m_prof = Float64[]
nom_e_prof = Float64[]; vac_e_prof = Float64[]
nom_pdg    = Int32[];   vac_pdg    = Int32[]


for (i, tree) in enumerate(trees)
    for row in tree

        E    = Float64(row.clusterEnergy)

        (E < energy_min || E > energy_max) && continue

        pdg    = Int32(row.pdgCode)
        ps     = Float64(row.profileStart)
        pd     = Float64(row.profileDiscrepancy)
        sso    = Float64(row.showerStartOffset)
        prms   = Float64(row.peakRms)
        rrms   = Float64(row.rmsRatio)
        #m_prof = Float64(row.measuredProfile)
        #e_prof = Float64(row.expectedProfile)

        if i==1
            push!(nom_ps, ps)
            push!(nom_pd, pd)
            push!(nom_e, E)
            push!(nom_sso, sso)
            push!(nom_prms, prms)
            push!(nom_rrms, rrms)
            push!(nom_pdg, pdg)
            #push!(nom_m_prof, m_prof)
            #push!(nom_e_prof, e_prof)
        elseif i==2
            push!(vac_ps, ps)
            push!(vac_pd, pd)
            push!(vac_e, E)
            push!(vac_sso, sso)
            push!(vac_prms, prms)
            push!(vac_rrms, rrms)
            push!(vac_pdg, pdg)
            #push!(vac_m_prof, m_prof)
            #push!(vac_e_prof, e_prof)
        else
            println("i=$(i), something has gone wrong")
        end
    end
end

# --- Do Stats/Calculations ---
println("Calculating stats")

num_photons_nom = length(filter(x->x==22, nom_pdg))
num_photons_vac = length(filter(x->x==22, vac_pdg))
frac_photons_nom = num_photons_nom/length(nom_pdg)
frac_photons_vac = num_photons_vac/length(vac_pdg) 
println("Fraction photons in nominal case: $(frac_photons_nom)")
println("Fraction photons in vacuum case: $(frac_photons_vac)")

frac_in_range_nom = (length(nom_e)/length(tree_nom))
frac_in_range_vac = (length(vac_e)/length(tree_vac))
println("Fraction in range $(energy_min) GeV to $(energy_max) in nominal case: $(frac_in_range_nom)")
println("Fraction in range $(energy_min) GeV to $(energy_max) in vacuum case: $(frac_in_range_vac)")

#=
println("Iterating over curves")
for (i, curve) in enumerate(nom_m_prof)
end 
=#

# --- Plot Things ---
const COLORS = Dict(
  :purple => "#101c4d", #purple
  :nominal => "#e38153", #orange
  :vacuum => "#d14382", #pink
  :blue => "#0052c5", #blue
  :gold => "#e3ab0f", #gold
  :yellow => "#ffe436", #yellow
  :text => "#eeeeee", #white
  :axis => "#eeeeee", #white
)

ALPHAS = Dict(
  :nominal => 0.6,
  :vacuum => 0.6,
)

function sep_plot(title, xlabel, nom_vals, vac_vals;
                  nbins=50, x_lo=nothing, x_hi=nothing, logy=false)
    all_vals = [nom_vals; vac_vals]
    xlo = isnothing(x_lo) ? minimum(all_vals) : x_lo
    xhi = isnothing(x_hi) ? maximum(all_vals) : x_hi
    edges = range(xlo, xhi, length=nbins+1)
 
    function norm_hist(vals)
        h = fit(Histogram, Float64.(vals), edges)
        s = sum(h.weights)
        s > 0 ? (h.weights ./ s) : h.weights
    end
 
    nom_w = norm_hist(nom_vals)
    vac_w = norm_hist(vac_vals)
    centers = collect(edges[1:end-1]) .+ step(edges) / 2
 
    fig = Figure(backgroundcolor = :transparent, size = (800, 600))
 
    ax_kwargs = (
        backgroundcolor  = :transparent,
        title            = title,
        xlabel           = xlabel,
        ylabel           = "Normalised entries / bin",
        titlecolor       = COLORS[:text],
        xlabelcolor      = COLORS[:text],
        ylabelcolor      = COLORS[:text],
        xticklabelcolor  = COLORS[:text],
        yticklabelcolor  = COLORS[:text],
        xgridvisible     = false,
        ygridvisible     = false,
        bottomspinecolor = COLORS[:axis],
        leftspinecolor   = COLORS[:axis],
        topspinecolor    = COLORS[:axis],
        rightspinecolor  = COLORS[:axis],
        xtickcolor       = COLORS[:axis],
        ytickcolor       = COLORS[:axis],
    )
 
    ax = logy ?
        Axis(fig[1, 1]; ax_kwargs..., yscale = log10) :
        Axis(fig[1, 1]; ax_kwargs...)
 
    # Overlaid, semi-transparent bars with solid outlines — matches the
    # original aesthetic (fillalpha=0.4, linewidth=1.7).
    barplot!(ax, centers, nom_w;
        width        = step(edges),
        color        = (COLORS[:nominal], ALPHAS[:nominal]),
        strokecolor  = COLORS[:nominal],
        strokewidth  = 1.7,
        gap          = 0,
        label        = "Nominal",
    )
    barplot!(ax, centers, vac_w;
        width        = step(edges),
        color        = (COLORS[:vacuum], ALPHAS[:vacuum]),
        strokecolor  = COLORS[:vacuum],
        strokewidth  = 1.7,
        gap          = 0,
        label        = "Vacuum",
    )
 
    if logy
        ylims!(ax, 1e-4, nothing)
    else
        ylims!(ax, 0, nothing)
    end
 
    axislegend(ax;
        position         = :rt,
        framevisible     = false,
        backgroundcolor  = :transparent,
        labelcolor       = COLORS[:text],
    )
 
    return fig
end


# --- Build all panels ---
all_e    = [nom_e;    vac_e]
all_pdg  = [nom_pdg;  vac_pdg]
all_pd   = [nom_pd;   vac_pd]
all_ps   = [nom_ps;   vac_ps]
all_sso  = [nom_sso;  vac_sso]
all_prms = [nom_prms; vac_prms]
all_rrms = [nom_rrms; vac_rrms]

# (filename_stem, Figure) — filename_stem becomes <output_dir>/<stem>.svg
pages = Tuple{String, Figure}[]
println("Adding Pages")

push!(pages, ("energy", sep_plot(
    "Energy: Nominal vs Vacuum Solenoid",
    "Energy [GeV]",
    nom_e, vac_e;
    nbins=50, x_lo=0.0, x_hi=min(maximum(all_e) * 1.05, 30),
)))

push!(pages, ("pdg", sep_plot(
    "PDG ID: Nominal vs Vacuum Solenoid",
    "ID",
    nom_pdg, vac_pdg;
    nbins=50, x_lo=0.0, x_hi=min(maximum(all_pdg) * 1.05, 30),
)))

push!(pages, ("profileDiscrepancy", sep_plot(
    "Profile Discrepancy: Nominal vs Vacuum Solenoid",
    L"\chi",
    nom_pd, vac_pd;
    nbins=50, x_lo=0.0, x_hi=min(maximum(all_pd) * 1.05, 30),
)))

push!(pages, ("profileStart", sep_plot(
    "Profile Start: Nominal vs Vacuum Solenoid",
    L"X_0",
    nom_ps, vac_ps;
    nbins=50, x_lo=0.0, x_hi=min(maximum(all_ps) * 1.05, 30),
)))

push!(pages, ("showerStartOffset", sep_plot(
    "Shower Start Offset (showerStartLayer - innerLayer): Nominal vs Vacuum Solenoid",
    "Pseudo Layer",
    nom_sso, vac_sso;
    nbins=50, x_lo=0.0, x_hi=min(maximum(all_sso) * 1.05, 30),
)))

push!(pages, ("peakRms", sep_plot(
    "Peak RMS: Nominal vs Vacuum Solenoid",
    "Peak RMS",
    nom_prms, vac_prms;
    nbins=50, x_lo=0.0, x_hi=min(maximum(all_prms) * 1.05, 30),
)))

push!(pages, ("rmsRatio", sep_plot(
    "RMS Ratio: Nominal vs Vacuum Solenoid",
    "RMS Ratio",
    nom_rrms, vac_rrms;
    nbins=50, x_lo=0.0, x_hi=min(maximum(all_rrms) * 1.05, 30),
)))

# --- Save ---
# SVG has no multi-page concept, so we write one file per variable.
mkpath(output_dir)
println("Saving $(length(pages)) SVG(s) to $output_dir/ ...")
for (stem, fig) in pages
    path = joinpath(output_dir, stem * ".svg")
    save(path, fig, backgroundcolor = :transparent)
    println("  wrote $path")
end
println("Done.")


