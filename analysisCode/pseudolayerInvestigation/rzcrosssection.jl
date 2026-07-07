#!/usr/bin/env julia
# rzCrossSection.jl  (v2)
#
# Produces, from the CSVs emitted by DumpPseudoLayerGeometry.cc /
# DumpSubDetectorEnvelopes.cc:
#   1. maia_rz_cross_section.pdf   - true-to-scale r-z, colored by pseudolayer
#   2. maia_ecal_x0_profile.pdf    - cumulative X0 vs pseudolayer (what the
#                                    PandoraPFA longitudinal profile integrates)
#   3. maia_material_budget.pdf    - OPTIONAL, from materialScan CSVs: the TRUE
#                                    radial material budget, where the coil shows up
#
# Deps:  ] add CSV DataFrames Plots

using CSV, DataFrames, Plots
gr()

# ---------------- configure your geometry cases here ----------------
cases = [
    (label="nominal", layers="maia_pseudolayers.csv", env="maia_envelopes.csv"),
    # (label="vacuum", layers="maia_pseudolayers_vac.csv", env="maia_envelopes_vac.csv"),
]

# Optional: true material-budget traces from `materialScan` (see header notes below).
# Schema: path_mm,cumX0,cumLambda
matscans = [
    # (label="nominal", file="matscan_nominal.csv"),
    # (label="vacuum",  file="matscan_vacuum.csv"),
]

# Solenoid extent for the annotation (NOT in the CSVs - it is a COIL DetElement).
# Set to `nothing` to hide. Replace with your real coil r/z.
solenoid = (rin=1750.0, rout=1857.0, halfz=2307.0)
# --------------------------------------------------------------------

load_case(c) = begin
    L = CSV.read(c.layers, DataFrame)
    E = CSV.read(c.env, DataFrame)
    (label=c.label, L=L, E=E, envd=Dict(String(r.subdetector) => r for r in eachrow(E)))
end
data = [load_case(c) for c in cases]

plcolor(pl, plmax) = cgrad(:viridis)[clamp(pl / plmax, 0.0, 1.0)]

# ---------------- Panel 1: r-z cross-section (first case) ----------------
function cross_section(d)
    L = d.L
    plmax = maximum(x -> x < 0 ? 0 : x, skipmissing(L.pseudolayer))
    plt = plot(size=(1150, 600), legend=false, aspect_ratio=:equal, framestyle=:box,
               xlabel="z [mm]", ylabel="r [mm]",
               title="MAIA r-z cross-section - $(d.label)  (color = pseudolayer)")
    if solenoid !== nothing
        s = solenoid
        plot!(plt, [-s.halfz, s.halfz, s.halfz, -s.halfz, -s.halfz],
                   [s.rin, s.rin, s.rout, s.rout, s.rin],
              seriestype=:shape, fillalpha=0.18, fillcolor=:red, linecolor=:red, lw=0.5)
        annotate!(plt, 0.0, (s.rin + s.rout) / 2, text("solenoid ~4 X0", 8, :red, :center))
    end
    for row in eachrow(L)
        (ismissing(row.pseudolayer) || row.pseudolayer < 0) && continue
        haskey(d.envd, String(row.subdetector)) || continue
        e = d.envd[String(row.subdetector)]
        c = plcolor(row.pseudolayer, plmax)
        if row.region == "barrel"
            plot!(plt, [-e.outerZ_mm, e.outerZ_mm], [row.position_mm, row.position_mm], color=c, lw=1.1)
        else
            plot!(plt, [ row.position_mm,  row.position_mm], [e.innerR_mm, e.outerR_mm], color=c, lw=1.1)
            plot!(plt, [-row.position_mm, -row.position_mm], [e.innerR_mm, e.outerR_mm], color=c, lw=1.1)
        end
    end
    plt
end

# ---------------- Panel 2: cumulative X0 vs pseudolayer (ECAL) ----------------
# This is the depth axis the LCShowerProfilePlugin actually integrates over.
# NOTE: for nominal vs vacuum these curves OVERLAP exactly - the coil X0 is not
# in the ECAL LayeredCalorimeterData. That overlap is the finding, not a bug.
function x0_profile()
    plt = plot(size=(950, 540), framestyle=:box, xlabel="pseudolayer",
               ylabel="cumulative X0", legend=:bottomright,
               title="ECAL longitudinal depth (calo layers only)")
    for d in data
        for (reg, ls) in (("ECAL_BARREL", :solid), ("ECAL_ENDCAP", :dash))
            sub = sort(filter(r -> String(r.subdetector) == reg &&
                                   !ismissing(r.pseudolayer) && r.pseudolayer >= 0, d.L),
                       :pseudolayer)
            nrow(sub) == 0 && continue
            tag = reg == "ECAL_BARREL" ? "barrel" : "endcap"
            plot!(plt, sub.pseudolayer, sub.cumulative_X0, ls=ls, lw=2, label="$(d.label) $tag")
        end
    end
    plt
end

# ---------------- Panel 3 (optional): true material budget ----------------
function material_budget()
    isempty(matscans) && return nothing
    plt = plot(size=(950, 540), framestyle=:box, xlabel="radial path length [mm]",
               ylabel="integrated X0", legend=:topleft,
               title="Material budget from materialScan (true, incl. coil)")
    if solenoid !== nothing
        vspan!(plt, [solenoid.rin, solenoid.rout], alpha=0.15, color=:red, label="solenoid")
    end
    for m in matscans
        df = CSV.read(m.file, DataFrame)
        plot!(plt, df.path_mm, df.cumX0, lw=2, label=m.label)
    end
    plt
end

savefig(cross_section(data[1]), "pseudolayerPlots/maia_rz_cross_section.pdf")
savefig(x0_profile(),           "pseudolayerPlots/maia_ecal_x0_profile.pdf")
mb = material_budget()
mb !== nothing && savefig(mb, "pseudolayerPlots/maia_material_budget.pdf")
println("wrote figures")
