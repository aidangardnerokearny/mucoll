#!/usr/bin/env julia
#=
plot_ecfit.jl — plots for the E_c fit (fit_plan.md, Deliverables), from the CSVs written by Ecfit.jl.
Usage (from this directory, after Ecfit.jl):  julia plot_ecfit.jl [fit_output_dir=.] [plot_dir=fit_output_dir] [ext=pdf]

  ecfit_scans.pdf         L_j(u) = <δ> vs E_c for each energy, normalised to its minimum, with the fitted minima
  ecfit_uhat_vs_E.pdf     fitted E_c per energy (±1σ bootstrap), constant-E_c and joint (α, p1) fits, 80 MeV default (F2)
  ecfit_mean_profiles.pdf <O_i>/E_ECAL per energy with P_i at the default and the adopted parameters, plus residuals
  ecfit_f3_classes.pdf    main-cluster fit by fragmentation class, Δu relative to the merged-shower fit (F3)
=#
using CSV, DataFrames, CairoMakie, Printf

dir = length(ARGS) >= 1 ? ARGS[1] : "."
pdir = length(ARGS) >= 2 ? ARGS[2] : dir
ext = length(ARGS) >= 3 ? ARGS[3] : "pdf"
mkpath(pdir)
res = CSV.read(joinpath(dir, "results.csv"), DataFrame; missingstring = "nan")
scans = CSV.read(joinpath(dir, "scans.csv"), DataFrame; missingstring = "nan")
mp = CSV.read(joinpath(dir, "mean_profiles.csv"), DataFrame; missingstring = "nan")

const INK = :black
set_theme!(Theme(fontsize = 13, textcolor = INK,
    Axis = (xgridcolor = (:black, 0.08), ygridcolor = (:black, 0.08), bottomspinecolor = INK, leftspinecolor = INK,
            topspinecolor = INK, rightspinecolor = INK, xtickcolor = INK, ytickcolor = INK,
            xticklabelcolor = INK, yticklabelcolor = INK, xlabelcolor = INK, ylabelcolor = INK, titlecolor = INK),
    Legend = (framecolor = (:black, 0.3), labelcolor = INK, titlecolor = INK)))

# One-hue sequential ramp (energy is ordered), mid→dark steps so every line is visible on white
const RAMP = ["#86b6ef", "#6da7ec", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#0d366b"]
const CAT = ["#2a78d6", "#eb6834", "#1baf7a"]           # categorical slots 1–3 (validated all-pairs)
const DEF_EC = 0.08

f1 = sort(res[res.fit .== "F1", :], :energy_gev)
labels = f1.name
energies = f1.energy_gev
elabel(e) = e < 1 ? @sprintf("%.1f GeV", e) : @sprintf("%g GeV", e)
comb = res[res.fit .== "F2_combined", :]
ecComb, σComb = isempty(comb) ? (NaN, NaN) : (comb.ec_gev[1], comb.sigma[1])

# ---------------- 1. scans ----------------
fig = Figure(size = (720, 470))
ax = Axis(fig[1, 1]; xscale = log10, xlabel = "E_c [GeV]", ylabel = "L(E_c) / L_min", yscale = log10,
          title = "Fixed-start discrepancy, merged photon showers (S1, vacuum solenoid)")
for (k, lab) in enumerate(labels)
    s = sort(scans[(scans.scope .== "sample") .& (scans.name .== lab), :], :u)
    lines!(ax, s.ec_gev, s.L ./ minimum(s.L); color = RAMP[k], linewidth = 2, label = elabel(energies[k]))
    scatter!(ax, [f1.ec_gev[k]], [1.0]; color = RAMP[k], markersize = 9, strokecolor = :white, strokewidth = 2)
end
vlines!(ax, [DEF_EC]; color = INK, linestyle = :dash, linewidth = 1.5, label = "default 80 MeV")
Legend(fig[1, 2], ax, "Photon energy")
save(joinpath(pdir, "ecfit_scans.$ext"), fig)

# ---------------- 2. û vs E ----------------
joint = res[res.fit .== "F2_joint", :]
αJ, p1J = isempty(joint) ? (NaN, NaN) : (joint.u_or_alpha[1], joint.ec_lo_or_p1[1])
adoptLine = filter(l -> startswith(l, "ADOPT"), readlines(joinpath(dir, "summary.txt")))
adoptedJoint = !isempty(adoptLine) && occursin("joint", adoptLine[1])
P0, P1DEF = 1.25, 0.5
fig = Figure(size = (680, 450))
mev = [5, 10, 20, 50, 100]
ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "Photon energy [GeV]", ylabel = "E_c [MeV]",
          title = "F1 per energy and F2 fits (S1, merged showers)",
          xticks = ([0.5, 1, 2, 5, 10, 25, 50, 100], ["0.5", "1", "2", "5", "10", "25", "50", "100"]),
          yticks = (mev, string.(mev)))
xr = [minimum(energies) / 1.5, maximum(energies) * 1.5]
if !isnan(ecComb)
    band!(ax, xr, fill(1000ecComb * exp(-σComb), 2), fill(1000ecComb * exp(σComb), 2); color = (CAT[1], 0.15))
    hlines!(ax, [1000ecComb]; color = CAT[1], linewidth = 2, label = @sprintf("constant E_c fit: %.1f MeV", 1000ecComb))
end
if !isnan(αJ)
    # E_c that reproduces the joint fit's a(E) = α + p1 ln E with the default p0, p1 (evaluated at E_clus = E)
    eg = exp.(range(log(xr[1]), log(xr[2]); length = 100))
    ecEq = eg .* exp.(-((αJ .+ p1J .* log.(eg)) .- P0) ./ P1DEF)
    lines!(ax, eg, 1000 .* ecEq; color = CAT[2], linewidth = 2,
           label = @sprintf("joint fit (p1 = %.3f) as E_c at p1 = 0.5%s", p1J, adoptedJoint ? " [adopted]" : ""))
end
hlines!(ax, [1000DEF_EC]; color = INK, linestyle = :dash, linewidth = 1.5, label = "default: 80 MeV")
errorbars!(ax, energies, 1000 .* f1.ec_gev, 1000 .* (f1.ec_gev .- f1.ec_lo_or_p1), 1000 .* (f1.ec_hi_or_sigma_p1 .- f1.ec_gev);
           color = INK, whiskerwidth = 8)
scatter!(ax, energies, 1000 .* f1.ec_gev; color = INK, markersize = 9, label = "F1 per energy (±1σ bootstrap)")
xlims!(ax, xr...); ylims!(ax, 7, 110)
axislegend(ax; position = :rb, labelsize = 11)
save(joinpath(pdir, "ecfit_uhat_vs_E.$ext"), fig)

# ---------------- 3. mean profiles ----------------
fig = Figure(size = (1100, 900))
ecFitLabel = adoptedJoint ? @sprintf("adopted joint fit (p1 = %.3f)", p1J) : @sprintf("adopted %.1f MeV", 1000ecComb)
for (k, lab) in enumerate(labels)
    r, c = fldmod1(k, 4)
    s = sort(mp[mp.sample .== lab, :], :bin)
    s = s[s.n .>= 0.5 * maximum(s.n), :]            # bins reached by at least half of the showers
    gl = fig[r, c] = GridLayout()
    a1 = Axis(gl[1, 1]; title = elabel(energies[k]), ylabel = k % 4 == 1 ? "<O_i> / E_ECAL" : "", xticklabelsvisible = false)
    a2 = Axis(gl[2, 1]; xlabel = "t [X0]", ylabel = k % 4 == 1 ? "data − P" : "")
    rowsize!(gl, 2, Relative(0.3)); rowgap!(gl, 4); linkxaxes!(a1, a2)
    scatter!(a1, s.t_upper_x0, s.mean_o; color = INK, markersize = 5, label = "data")
    lines!(a1, s.t_upper_x0, s.mean_p_default; color = CAT[2], linewidth = 2, label = "P, default 80 MeV")
    lines!(a1, s.t_upper_x0, s.mean_p_fit; color = CAT[1], linewidth = 2, label = "P, " * ecFitLabel)
    hlines!(a2, [0.0]; color = (:black, 0.4), linewidth = 1)
    lines!(a2, s.t_upper_x0, s.mean_o .- s.mean_p_default; color = CAT[2], linewidth = 2)
    lines!(a2, s.t_upper_x0, s.mean_o .- s.mean_p_fit; color = CAT[1], linewidth = 2)
    k == 1 && (global legAx = a1)
end
Legend(fig[3, :], legAx; orientation = :horizontal, framevisible = false)
Label(fig[0, :], "Mean longitudinal profile of merged photon showers vs expected profile (S1)"; fontsize = 15, color = INK)
save(joinpath(pdir, "ecfit_mean_profiles.$ext"), fig)

# ---------------- 4. F3 by fragmentation class ----------------
fig = Figure(size = (640, 430))
ax = Axis(fig[1, 1]; xscale = log10, xlabel = "Photon energy [GeV]", ylabel = "û(main cluster) − û(merged)",
          title = "F3: fragmentation bias by class (A < 1%, B 1–5%, C > 5%)")
hlines!(ax, [0.0]; color = (:black, 0.4), linewidth = 1)
offsets = Dict('A' => 0.94, 'B' => 1.0, 'C' => 1.06)
for (j, cls) in enumerate(('A', 'B', 'C'))
    d = res[res.fit .== "F3_$cls", :]
    isempty(d) && continue
    uref = [f1.u_or_alpha[findfirst(==(n), f1.name)] for n in d.name]
    x = d.energy_gev .* offsets[cls]
    errorbars!(ax, x, d.u_or_alpha .- uref, d.sigma; color = CAT[j], whiskerwidth = 6)
    scatterlines!(ax, x, d.u_or_alpha .- uref; color = CAT[j], markersize = 9, linewidth = 1.5,
                  label = "class $cls")
end
axislegend(ax; position = :lb)
save(joinpath(pdir, "ecfit_f3_classes.$ext"), fig)

println("wrote ecfit_{scans,uhat_vs_E,mean_profiles,f3_classes}.$ext in ", abspath(pdir))
