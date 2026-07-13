#!/usr/bin/env julia
# rzCrossSectionCustom.jl
#
# CairoMakie rewrite of rzcrosssection.jl, styled to match profileCompCustom.jl
# (transparent background, custom font/color scheme, .svg output).
#
# Produces, from the CSVs emitted by DumpPseudoLayerGeometry.cc /
# DumpSubDetectorEnvelopes.cc:
#   1. maia_rz_cross_section.svg - true-to-scale r-z, colored by pseudolayer
#   2. maia_ecal_x0_profile.svg  - cumulative X0 vs pseudolayer (the depth axis
#                                  the LCShowerProfilePlugin actually integrates over)
#
# Deps:  ] add CSV DataFrames CairoMakie

using CSV
using DataFrames
using CairoMakie

set_theme!(fonts = (
    regular     = "/home/aidangardnerokearny/.local/share/fonts/OpenSans/OpenSans-Regular.ttf",
    bold        = "/home/aidangardnerokearny/.local/share/fonts/OpenSans/OpenSans-Bold.ttf",
    italic      = "/home/aidangardnerokearny/.local/share/fonts/OpenSans/OpenSans-Italic.ttf",
    bold_italic = "/home/aidangardnerokearny/.local/share/fonts/OpenSans/OpenSans-BoldItalic.ttf",
))

# --- Arguments ---
layers_nom = length(ARGS) >= 1 ? ARGS[1] : "maia_pseudolayers.csv"
env_nom    = length(ARGS) >= 2 ? ARGS[2] : "maia_envelopes.csv"
output_dir = length(ARGS) >= 3 ? ARGS[3] : "pseudolayerPlots"
layers_vac = length(ARGS) >= 4 ? ARGS[4] : nothing
env_vac    = length(ARGS) >= 5 ? ARGS[5] : nothing

# Solenoid extent for the annotation (NOT in the CSVs - it is a COIL DetElement).
# Set to `nothing` to hide. Replace with your real coil r/z.
solenoid = (rin=1750.0, rout=1857.0, halfz=2307.0)

# --- Load Cases ---
function load_case(label, layers_file, env_file)
    println("Loading $layers_file / $env_file ...")
    L = CSV.read(layers_file, DataFrame)
    E = CSV.read(env_file, DataFrame)
    envd = Dict(String(r.subdetector) => r for r in eachrow(E))
    (label=label, L=L, E=E, envd=envd)
end

cases = [load_case("nominal", layers_nom, env_nom)]
if layers_vac !== nothing && env_vac !== nothing
    push!(cases, load_case("vacuum", layers_vac, env_vac))
end

# --- Plot Style ---
const COLORS = Dict(
  :purple  => "#101c4d", #purple
  :nominal => "#e38153", #orange
  :vacuum  => "#d14382", #pink
  :blue    => "#0052c5", #blue
  :gold    => "#e3ab0f", #gold
  :yellow  => "#ffe436", #yellow
  :text    => "#eeeeee", #white
  :axis    => "#eeeeee", #white
)

const CASE_COLOR = Dict("nominal" => COLORS[:nominal], "vacuum" => COLORS[:vacuum])

function styled_axis(fig_pos; xlabel, ylabel, kwargs...)
    Axis(fig_pos;
        backgroundcolor    = :transparent,
        xlabel             = xlabel,
        ylabel             = ylabel,
        xlabelcolor        = COLORS[:text],
        ylabelcolor        = COLORS[:text],
        xticklabelcolor    = COLORS[:text],
        yticklabelcolor    = COLORS[:text],
        xgridvisible       = false,
        ygridvisible       = false,
        bottomspinecolor   = COLORS[:axis],
        leftspinecolor     = COLORS[:axis],
        topspinecolor      = COLORS[:axis],
        rightspinecolor    = COLORS[:axis],
        xtickcolor         = COLORS[:axis],
        ytickcolor         = COLORS[:axis],
        xminorticksvisible = true,
        yminorticksvisible = true,
        xminortickcolor    = COLORS[:axis],
        yminortickcolor    = COLORS[:axis],
        xminorticks        = IntervalsBetween(5),
        yminorticks        = IntervalsBetween(5),
        kwargs...,
    )
end

# --- Panel 1: r-z cross-section (true to scale, colored by pseudolayer) ---
function cross_section(d)
    L = d.L
    plmax = maximum(x -> (ismissing(x) || x < 0) ? 0 : x, L.pseudolayer)

    fig = Figure(backgroundcolor = :transparent, size = (1150, 650))
    ax = styled_axis(fig[1, 1]; xlabel = "z [mm]", ylabel = "r [mm]", aspect = DataAspect())

    if solenoid !== nothing
        s = solenoid
        poly!(ax, Point2f[(-s.halfz, s.rin), (s.halfz, s.rin), (s.halfz, s.rout), (-s.halfz, s.rout)];
              color = (:red, 0.18), strokecolor = :red, strokewidth = 0.5)
        text!(ax, 0.0, (s.rin + s.rout) / 2;
              text = "solenoid ~4 X0", color = :red, fontsize = 12, align = (:center, :center))
    end

    for row in eachrow(L)
        (ismissing(row.pseudolayer) || row.pseudolayer < 0) && continue
        haskey(d.envd, String(row.subdetector)) || continue
        e = d.envd[String(row.subdetector)]
        pl = Float64(row.pseudolayer)
        if row.region == "barrel"
            lines!(ax, [-e.outerZ_mm, e.outerZ_mm], [row.position_mm, row.position_mm];
                   color = pl, colormap = :viridis, colorrange = (0, plmax), linewidth = 1.4)
        else
            lines!(ax, [row.position_mm, row.position_mm], [e.innerR_mm, e.outerR_mm];
                   color = pl, colormap = :viridis, colorrange = (0, plmax), linewidth = 1.4)
            lines!(ax, [-row.position_mm, -row.position_mm], [e.innerR_mm, e.outerR_mm];
                   color = pl, colormap = :viridis, colorrange = (0, plmax), linewidth = 1.4)
        end
    end

    Colorbar(fig[1, 2]; colormap = :viridis, limits = (0, plmax), label = "pseudolayer",
             labelcolor = COLORS[:text], ticklabelcolor = COLORS[:text],
             tickcolor = COLORS[:axis], spinewidth = 0)

    Label(fig[0, :], "MAIA r-z cross-section — $(d.label)"; color = COLORS[:text], fontsize = 20)

    fig
end

# --- Panel 2: cumulative X0 vs pseudolayer (ECAL longitudinal depth) ---
# NOTE: for nominal vs vacuum these curves OVERLAP exactly - the coil X0 is not
# in the ECAL LayeredCalorimeterData. That overlap is the finding, not a bug.
function x0_profile(cases)
    fig = Figure(backgroundcolor = :transparent, size = (950, 560))
    ax = styled_axis(fig[1, 1]; xlabel = "pseudolayer", ylabel = "cumulative X0")

    for d in cases
        col = get(CASE_COLOR, d.label, COLORS[:blue])
        for (reg, ls, tag) in (("ECAL_BARREL", :solid, "barrel"), ("ECAL_ENDCAP", :dash, "endcap"))
            sub = sort(filter(r -> String(r.subdetector) == reg &&
                                   !ismissing(r.pseudolayer) && r.pseudolayer >= 0, d.L),
                       :pseudolayer)
            nrow(sub) == 0 && continue
            lines!(ax, sub.pseudolayer, sub.cumulative_X0;
                   color = col, linestyle = ls, linewidth = 2, label = "$(d.label) $tag")
        end
    end

    axislegend(ax; position = :rb, framevisible = false, backgroundcolor = :transparent,
               labelcolor = COLORS[:text])
    Label(fig[0, :], "ECAL longitudinal depth (calo layers only)"; color = COLORS[:text], fontsize = 20)
    fig
end

# --- Save ---
mkpath(output_dir)
pages = Tuple{String, Figure}[]
push!(pages, ("maia_rz_cross_section", cross_section(cases[1])))
push!(pages, ("maia_ecal_x0_profile",  x0_profile(cases)))

println("Saving $(length(pages)) SVG(s) to $output_dir ...")
for (stem, fig) in pages
    path = joinpath(output_dir, stem * ".svg")
    save(path, fig, backgroundcolor = :transparent)
    println("  wrote $path")
end
println("Done.")
