"""
plot_shower_profiles.jl
 
Plots measured and expected longitudinal shower energy-deposition curves
from shower_profiles.root using Julia (UnROOT + Plots).
 
Pages in output PDF:
  1..N  Individual cluster measured vs expected profiles (one per cluster)
  N+1   Overlay of selected clusters across energy range
  N+2   Average measured vs average expected profile
  N+3   profileStart distribution
  N+4   profileDiscrepancy distribution
  N+5   Pseudo-layer distributions (innerLayer vs showerStartLayer)
  N+6   2D cluster energy vs nCaloHits
  N+7   1D cluster energy distribution (log-y)
 
Usage
-----
    julia plot_shower_profiles.jl [input.root] [output.pdf] [max_clusters] [e_min] [e_max]
 
    Defaults: shower_profiles.root  shower_profiles.pdf  6  0.5  1e9
"""
 
using UnROOT
using Plots
using StatsBase
using Printf
using Statistics
 
# ── Argument parsing ───────────────────────────────────────────────────────────
input_file   = length(ARGS) >= 1 ? ARGS[1] : "shower_profiles.root"
output_file  = length(ARGS) >= 2 ? ARGS[2] : "shower_profiles.pdf"
max_clusters = length(ARGS) >= 3 ? parse(Int,     ARGS[3]) : 6
energy_min   = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.5
energy_max   = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1e9
 
const BIN_WIDTH = 0.5f0   # radiation lengths per bin
 
# ── Global plot defaults ───────────────────────────────────────────────────────
default(
    grid        = true,
    gridalpha   = 0.3,
    framestyle  = :box,
    legend      = :topright,
    legendframe = false,   # no legend box
    fontfamily  = "sans-serif",
    titlefont   = 10,
    labelfont   = 9,
    tickfont    = 8,
)
 
# ── Open ROOT file ─────────────────────────────────────────────────────────────
println("Opening $input_file ...")
f    = ROOTFile(input_file)
tree = f["ShowerProfiles"]
println("  $(length(tree)) entries in 'ShowerProfiles'")
 
# ── Read branches ──────────────────────────────────────────────────────────────
cluster_energy      = LazyTree(f, "ShowerProfiles", ["clusterEnergy"])        |> x -> x.clusterEnergy
n_calo_hits         = LazyTree(f, "ShowerProfiles", ["nCaloHits"])            |> x -> x.nCaloHits
profile_start       = LazyTree(f, "ShowerProfiles", ["profileStart"])         |> x -> x.profileStart
profile_disc        = LazyTree(f, "ShowerProfiles", ["profileDiscrepancy"])   |> x -> x.profileDiscrepancy
n_profile_bins      = LazyTree(f, "ShowerProfiles", ["nProfileBins"])         |> x -> x.nProfileBins
inner_layer         = LazyTree(f, "ShowerProfiles", ["innerLayer"])           |> x -> x.innerLayer
shower_start_layer  = LazyTree(f, "ShowerProfiles", ["showerStartLayer"])     |> x -> x.showerStartLayer
measured_profile    = LazyTree(f, "ShowerProfiles", ["measuredProfile"])      |> x -> x.measuredProfile
expected_profile    = LazyTree(f, "ShowerProfiles", ["expectedProfile"])      |> x -> x.expectedProfile
 
n_entries = length(cluster_energy)
 
# ── Collect data ───────────────────────────────────────────────────────────────
struct Cluster
    energy       :: Float32
    profile_start:: Float32
    profile_disc :: Float32
    inner_layer  :: UInt32
    shower_start :: UInt32
    n_bins       :: Int32
    n_hits       :: UInt32
    measured     :: Vector{Float32}
    expected     :: Vector{Float32}
end
 
clusters     = Cluster[]
all_energies = Float32[]
all_nhits    = UInt32[]
 
for i in 1:n_entries
    E    = cluster_energy[i]
    nhit = n_calo_hits[i]
    push!(all_energies, E)
    push!(all_nhits,    nhit)
 
    nb = n_profile_bins[i]
    nb == 0 && continue
    ps = profile_start[i]
    ps <= 0 && continue
    (E < energy_min || E > energy_max) && continue
 
    meas = collect(Float32, measured_profile[i])
    expt = collect(Float32, expected_profile[i])
 
    push!(clusters, Cluster(
        E, ps, profile_disc[i],
        inner_layer[i], shower_start_layer[i],
        nb, nhit, meas, expt
    ))
end
 
println("$(length(all_energies)) total clusters in tree")
println("$(length(clusters)) clusters pass selection cuts")
isempty(clusters) && error("No clusters pass the selection cuts.")
 
# ── Helper: t-axis for a cluster ──────────────────────────────────────────────
t_axis(cl) = Float32.(0:cl.n_bins-1) .* BIN_WIDTH
 
# ── Collect all plots into a vector for PDF output ────────────────────────────
pages = Plots.Plot[]
 
# ── Per-cluster pages ─────────────────────────────────────────────────────────
println("Building $(length(clusters)) individual cluster pages...")
for (idx, cl) in enumerate(clusters)
    t = t_axis(cl)
    p = plot(t, cl.measured;
        label     = "Measured",
        color     = :royalblue,
        lw        = 2,
        linestyle = :solid,
        xlabel    = "Depth (radiation lengths X₀)",
        ylabel    = "Energy deposit [GeV / bin]",
        title     = @sprintf("Cluster %d  E=%.3f GeV  profileStart=%.2f X₀  χ=%.4f",
                             idx, cl.energy, cl.profile_start, cl.profile_disc),
    )
    plot!(p, t, cl.expected;
        label     = "Expected",
        color     = :crimson,
        lw        = 2,
        linestyle = :dash,
    )
    push!(pages, p)
end
 
# ── Overlay of selected clusters ──────────────────────────────────────────────
clusters_sorted = sort(clusters, by = cl -> cl.energy)
n_show   = min(max_clusters, length(clusters_sorted))
step     = max(1, length(clusters_sorted) ÷ n_show)
selected = clusters_sorted[1:step:end][1:n_show]
 
palette_colors = [:royalblue, :crimson, :forestgreen,
                  :darkorange, :mediumpurple, :teal,
                  :hotpink, :saddlebrown]
 
p_overlay = plot(;
    xlabel = "Depth (radiation lengths X₀)",
    ylabel = "Energy deposit [GeV / bin]",
    title  = "Longitudinal shower profiles: measured (solid) vs expected (dashed)",
)
for (idx, cl) in enumerate(selected)
    t   = t_axis(cl)
    col = palette_colors[mod1(idx, length(palette_colors))]
    lbl = @sprintf("E=%.2f GeV  χ=%.3f", cl.energy, cl.profile_disc)
    plot!(p_overlay, t, cl.measured; label = lbl,  color = col, lw = 2, linestyle = :solid)
    plot!(p_overlay, t, cl.expected; label = false, color = col, lw = 2, linestyle = :dash)
end
push!(pages, p_overlay)
 
# ── Average measured vs average expected ──────────────────────────────────────
max_bins  = maximum(cl.n_bins for cl in clusters)
sum_meas  = zeros(Float64, max_bins)
sum_exp   = zeros(Float64, max_bins)
bin_count = zeros(Int,     max_bins)
 
for cl in clusters
    for i in 1:cl.n_bins
        sum_meas[i]  += cl.measured[i]
        sum_exp[i]   += cl.expected[i]
        bin_count[i] += 1
    end
end
 
valid    = bin_count .> 0
t_avg    = Float32.(0:max_bins-1)[valid] .* BIN_WIDTH
avg_meas = (sum_meas ./ max.(bin_count, 1))[valid]
avg_exp  = (sum_exp  ./ max.(bin_count, 1))[valid]
 
p_avg = plot(t_avg, avg_meas;
    label     = "Mean measured",
    color     = :royalblue,
    lw        = 3,
    linestyle = :solid,
    xlabel    = "Depth (radiation lengths X₀)",
    ylabel    = "Mean energy deposit [GeV / bin]",
    title     = "Average longitudinal profile ($(length(clusters)) clusters)",
)
plot!(p_avg, t_avg, avg_exp;
    label     = "Mean expected",
    color     = :crimson,
    lw        = 3,
    linestyle = :dash,
)
push!(pages, p_avg)
 
# ── profileStart distribution ──────────────────────────────────────────────────
ps_vals = [cl.profile_start for cl in clusters]
p_ps = histogram(ps_vals;
    bins    = 40,
    color   = :steelblue,
    alpha   = 0.6,
    label   = false,
    xlabel  = "profileStart (X₀)",
    ylabel  = "Clusters / bin",
    title   = "Shower profile start distribution",
)
push!(pages, p_ps)
 
# ── profileDiscrepancy distribution ───────────────────────────────────────────
pd_vals = [cl.profile_disc for cl in clusters]
p_pd = histogram(pd_vals;
    bins    = 40,
    color   = :tomato,
    alpha   = 0.6,
    label   = false,
    xlabel  = "profileDiscrepancy (normalised)",
    ylabel  = "Clusters / bin",
    title   = "Profile discrepancy distribution",
)
push!(pages, p_pd)
 
# ── Pseudo-layer distributions ────────────────────────────────────────────────
il_vals  = Float64[cl.inner_layer  for cl in clusters]
ssl_vals = Float64[cl.shower_start for cl in clusters]
layer_max = maximum([il_vals; ssl_vals]) * 1.05
 
p_layers = histogram(il_vals;
    bins    = range(0, layer_max, length = 80),
    color   = :royalblue,
    alpha   = 0.4,
    label   = "Inner pseudo-layer (cluster start)",
    xlabel  = "Pseudo-layer",
    ylabel  = "Clusters / bin",
    title   = "Cluster and shower start pseudo-layers",
)
histogram!(p_layers, ssl_vals;
    bins    = range(0, layer_max, length = 80),
    color   = :crimson,
    alpha   = 0.4,
    label   = "Shower start pseudo-layer",
)
push!(pages, p_layers)
 
# ── 2D cluster energy vs nCaloHits ────────────────────────────────────────────
e_max    = maximum(all_energies) * 1.05
hits_max = maximum(all_nhits)    * 1.05
ref_slope = 30.0
 
p_2d = histogram2d(Float64.(all_energies), Float64.(all_nhits);
    bins      = (100, 100),
    xlabel    = "Cluster EM energy [GeV]",
    ylabel    = "Number of calo hits",
    title     = "Cluster energy vs nCaloHits",
    color     = :viridis,
    colorbar  = true,
    label     = false,
)
plot!(p_2d, [0.0, e_max], [0.0, ref_slope * e_max];
    label     = "Reference: $(Int(ref_slope)) hits/GeV",
    color     = :red,
    lw        = 2,
    linestyle = :dash,
)
push!(pages, p_2d)
 
# ── 1D cluster energy distribution (log-y) ────────────────────────────────────
p_edist = histogram(Float64.(all_energies);
    bins    = 100,
    color   = :steelblue,
    alpha   = 0.6,
    label   = false,
    xlabel  = "Cluster EM energy [GeV]",
    ylabel  = "Clusters / bin",
    title   = "Cluster energy distribution (all clusters)",
    yscale  = :log10,
    ylims   = (0.5, :auto),
)
vline!(p_edist, [250.0, 1000.0];
    label     = "250-1000 GeV gun range",
    color     = :red,
    lw        = 2,
    linestyle = :dash,
)
push!(pages, p_edist)
 
# ── Save all pages to PDF ──────────────────────────────────────────────────────
println("Saving $(length(pages)) pages to $output_file ...")
savefig(pages[1], output_file)   # creates/overwrites the file
for p in pages[2:end]
    savefig(p, output_file)      # Plots.jl appends when backend is PDF
end
 
# ── Summary ───────────────────────────────────────────────────────────────────
energies  = [cl.energy        for cl in clusters]
ps_list   = [cl.profile_start for cl in clusters]
pd_list   = [cl.profile_disc  for cl in clusters]
n_sub_gev = count(e -> e < 1.0, all_energies)
 
println("\nSummary statistics ($(length(clusters)) selected clusters):")
println("  Total clusters in tree : $(length(all_energies))")
println("  Sub-GeV clusters       : $n_sub_gev ($(round(100*n_sub_gev/length(all_energies), digits=1))%)")
println("  Energy range (selected): $(round(minimum(energies), digits=2)) - $(round(maximum(energies), digits=2)) GeV")
println("  Mean profileStart      : $(round(mean(ps_list), digits=3)) ± $(round(std(ps_list), digits=3)) X₀")
println("  Mean discrepancy       : $(round(mean(pd_list), digits=4)) ± $(round(std(pd_list), digits=4))")
println("Saved $output_file")
 

