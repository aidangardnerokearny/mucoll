#=
averageProfileCompCustom.jl

Overlays the mean measured longitudinal shower profile for three geometries:
  - nominal
  - vacuum solenoid
  - all vacuum (everything before the ECAL is vacuum)

Averaging follows averageProfileCustom.py: each cluster's measured profile is
normalised to its own sum, then averaged bin-by-bin over the clusters that
have that bin.

Usage
-----
    julia averageProfileCompCustom.jl [nominal.root] [vacuum.root] [all_vacuum.root]
                                      [output_dir] [energy_min] [energy_max]
=#

using UnROOT
using CairoMakie

# --- Arguments ---
output_base = "/home/aidangardnerokearny/mucoll/jobs/output"
file_nom    = length(ARGS) >= 1 ? ARGS[1] : joinpath(output_base, "shower_profiles_nominal.root")
file_vac    = length(ARGS) >= 2 ? ARGS[2] : joinpath(output_base, "shower_profiles_vacuum.root")
file_allvac = length(ARGS) >= 3 ? ARGS[3] : joinpath(output_base, "up_vac", "shower_profiles.root")
output_dir  = length(ARGS) >= 4 ? ARGS[4] : "customPDFs"
energy_min  = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 0.5
energy_max  = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 1e9

const BIN_WIDTH = 0.5   # radiation lengths per bin (must match PandoraSettings)

# --- Average measured profile for one file ---
function average_profile(file)
    println("Opening $file...")
    tree = LazyTree(file, "ShowerProfiles")
    println("  $(length(tree)) entries")

    sum_meas = Float64[]
    counts   = Int[]
    n_pass   = 0

    for row in tree
        n_bins = Int(row.nProfileBins)
        n_bins == 0 && continue
        row.profileStart <= 0 && continue
        E = Float64(row.clusterEnergy)
        (E < energy_min || E > energy_max) && continue

        meas  = Float64.(row.measuredProfile[1:n_bins])
        total = sum(meas)
        total = total == 0 ? 1.0 : total

        if n_bins > length(sum_meas)
            append!(sum_meas, zeros(n_bins - length(sum_meas)))
            append!(counts,   zeros(Int, n_bins - length(counts)))
        end
        sum_meas[1:n_bins] .+= meas ./ total
        counts[1:n_bins]   .+= 1
        n_pass += 1
    end

    println("  $n_pass clusters pass selection cuts")
    n_pass == 0 && error("No clusters pass the selection cuts in $file")

    keep = counts .> 0
    t    = ((0:length(counts)-1) .* BIN_WIDTH)[keep]
    avg  = (sum_meas ./ max.(counts, 1))[keep]
    return t, avg
end

t_nom,    avg_nom    = average_profile(file_nom)
t_vac,    avg_vac    = average_profile(file_vac)
t_allvac, avg_allvac = average_profile(file_allvac)

# --- Plot ---
fig = Figure()
ax  = Axis(fig[1, 1]; xlabel = "Depth (X0)", ylabel = "Mean fractional energy deposit / bin")

lines!(ax, t_nom,    avg_nom;    label = "Nominal")
lines!(ax, t_vac,    avg_vac;    label = "Vacuum solenoid")
lines!(ax, t_allvac, avg_allvac; label = "All vacuum")
axislegend(ax)

# --- Save ---
mkpath(output_dir)
path = joinpath(output_dir, "average_profile_comparison.pdf")
save(path, fig)
println("Wrote $path")
