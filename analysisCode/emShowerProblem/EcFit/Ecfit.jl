#!/usr/bin/env julia
#=
ecFit.jl — fit LCShowerProfilePlugin's LongProfileCriticalEnergy (E_c) by minimising the
fixed-start profile discrepancy δ, following the E_c fit procedure note (ec_fit_procedure.pdf).

Input:  ROOT files written by ShowerProfileMonitoringAlgorithm (trees ShowerProfiles, Hits, Truth, Config).
Usage:  julia -t auto --project=. ecFit.jl ecfit.toml   (first time: julia --project=. -e "using Pkg; Pkg.instantiate()")
Needs:  UnROOT, SpecialFunctions (everything else is stdlib).

What it does, in order (section numbers refer to the procedure note):
  1. Reads every sample, rebuilds each photon's merged-hit-set profile with the plugin's exact rules (§2).
  2. Closure test (§7.2): rebuilds every dumped cluster from its own hits and requires agreement with
     the production profileStart / profileDiscrepancy / δ. Aborts if the pass fraction is too low.
  3. Selection (§5): truth endpoint, fiducial |cosθ|, valid profile, containment percentile.
  4. F1: per-sample fit of u = ln E_c minimising L(u) = <δ(u)> over merged showers, bootstrap σ.
  5. F2: per group, constancy test of û vs energy, combined equal-weight fit (adopted if constant),
     inverse-variance combination (systematic variant), and the joint (α, p1) fit (adopted if not).
  6. F3: main-cluster fit with production inputs, split by fragmentation class.
  7. F4: χ² fit of the mean normalised profile (L1-vs-mean comparison).
  F5 (electrons) and F6 (angles) are just extra samples/groups: give them their own `group`.

Conventions: u = ln(E_c / GeV). All energies in GeV, lengths in mm, depths in X0.
=#

using UnROOT
using SpecialFunctions: logabsgamma, gamma_inc
using Statistics
using Random
using Printf
using TOML

# =====================================================================================================================
# Configuration
# =====================================================================================================================

"""LCShowerProfilePlugin settings as dumped in the Config tree (types match the C++ members)."""
struct PluginConfig
    binWidth::Float32
    nBins::Int
    minCosAngle::Float32
    criticalEnergy::Float32
    parameter0::Float32
    parameter1::Float32
    maxDifference::Float32
    pseudoLayerAtIp::Int
end

const DEFAULTS = Dict(
    "output_dir" => "ecfit_out",
    "primary_group" => "S1",
    "selection" => Dict(
        "gun_pdg" => 22,
        "require_single_root" => true,
        "require_endpoint" => true,
        "min_conversion_radius_mm" => 0.0,   # barrel: truth endpoint radius must exceed this
        "max_abs_costheta" => 1.0,           # barrel fiducial
        "containment_percentile" => 1.0,     # drop events below this percentile of E_reco/E_true
        "gun_weight_min" => 0.5,             # hit belongs to the photon if gunWeightFraction >= this
    ),
    "fit" => Dict(
        "direction" => "truth",              # "truth" or "mainfit" (main cluster's fit direction)
        "ec_min_gev" => 1.0e-3,
        "ec_max_gev" => 1.0,
        "n_coarse" => 60,
        "fine_halfwidth" => 0.75,            # half-width in u of the bootstrap grid around û
        "fine_step" => 0.005,                # u step of the bootstrap grid (parabolic interpolation below this)
        "n_bootstrap" => 500,
        "n_bootstrap_joint" => 50,
        "seed" => 20261115,
        "fragmentation_edges" => [0.01, 0.05],
        "min_events_per_class" => 200,
        "f4_min_bin_occupancy" => 0.9,       # F4 uses bins populated in at least this fraction of showers
        "constancy_pvalue" => 0.05,
        "constancy_min_spread_u" => 0.02,    # a trend must also exceed this spread in ln Ec (2%) to count
    ),
    "closure" => Dict(
        "enabled" => true,
        "rel_tol" => 1.0e-5,
        "abs_tol_gev" => 1.0e-6,
        "require_pass_fraction" => 1.0,
        "max_clusters_per_sample" => 0,      # 0 = all
    ),
)

function merge_defaults(cfg::Dict, defaults::Dict)
    out = deepcopy(defaults)
    for (k, v) in cfg
        if haskey(out, k) && out[k] isa Dict && v isa Dict
            out[k] = merge_defaults(v, out[k])
        else
            out[k] = v
        end
    end
    return out
end

"""Expand a single '*' / '?' glob in the file-name part of a path (no extra dependency)."""
function expand_glob(pattern::AbstractString)
    dir, base = dirname(pattern), basename(pattern)
    occursin(r"[*?]", base) || return [pattern]
    rx = Regex("^" * replace(base, "." => "\\.", "*" => ".*", "?" => ".") * "\$")
    d = isempty(dir) ? "." : dir
    return sort([joinpath(dir, f) for f in readdir(d) if occursin(rx, f)])
end

# =====================================================================================================================
# Plugin replica: profile builder, expected profile, offset scan (LCShowerProfilePlugin::CalculateLongitudinalProfile)
# =====================================================================================================================

"""Columns of the Hits tree for one file."""
struct HitCols
    event::Vector{UInt32}
    cluster::Vector{Int32}
    isolated::Vector{Bool}
    pseudoLayer::Vector{UInt32}
    hitType::Vector{Int32}
    fine::Vector{Bool}
    em::Vector{Float32}
    x0::Vector{Float32}
    nx::Vector{Float32}
    ny::Vector{Float32}
    nz::Vector{Float32}
    w::Vector{Float32}
end

"""CartesianVector::GetCosOpeningAngle in Float32; returns NaN where Pandora would throw."""
function cos_opening(ax::Float32, ay::Float32, az::Float32, d::NTuple{3,Float32})
    magsq = (ax * ax + ay * ay + az * az) * (d[1] * d[1] + d[2] * d[2] + d[3] * d[3])
    magsq < eps(Float32) && return NaN32
    c = (ax * d[1] + ay * d[2] + az * d[3]) / sqrt(magsq)
    return clamp(c, -1.0f0, 1.0f0)
end

"""
    build_profile(h, idx, dir, cfg) -> (profile, eCalEnergy, profileEndBin, ok)

Step 1 of CalculateLongitudinalProfile for the hits `idx` (ordered, i.e. non-isolated, hits of one cluster
or of a merged set) and shower direction `dir`. `profile` has cfg.nBins bins; only 1:profileEndBin are used.
ok = false where the plugin returns early (coarse inner layer) or throws.
"""
function build_profile(h::HitCols, idx::AbstractVector{<:Integer}, dir::NTuple{3,Float32}, cfg::PluginConfig)
    profile = zeros(Float32, cfg.nBins)
    isempty(idx) && return (profile, 0.0f0, 0, false)

    ord = sort(collect(idx); by = i -> h.pseudoLayer[i], alg = Base.Sort.MergeSort)
    inner, outer = Int(h.pseudoLayer[ord[1]]), Int(h.pseudoLayer[ord[end]])

    # GetInnerLayerHitType: the hit type with the largest summed energy in the inner layer.
    # (Pandora uses hadronic energy; EM energy is the proxy here — identical for single-type layers.)
    typeE = Dict{Int32,Float32}()
    typeFine = Dict{Int32,Bool}()
    for i in ord
        h.pseudoLayer[i] == inner || break
        typeE[h.hitType[i]] = get(typeE, h.hitType[i], 0.0f0) + h.em[i]
        typeFine[h.hitType[i]] = h.fine[i]
    end
    innerType = first(sort(collect(typeE); by = last, rev = true))[1]
    typeFine[innerType] || return (profile, 0.0f0, 0, false)

    bw = cfg.binWidth
    eCal, nX0, nX0Last = 0.0f0, 0.0f0, 0.0f0
    k = 1
    for iLayer in inner:outer
        # hits of this layer are ord[k:kEnd]
        kEnd = k - 1
        while kEnd < length(ord) && h.pseudoLayer[ord[kEnd+1]] == iLayer
            kEnd += 1
        end
        if kEnd < k
            nX0 += nX0Last
            continue
        end
        nHits = kEnd - k + 1

        isFine, eLayer, x0Layer = true, 0.0f0, 0.0f0
        for j in k:kEnd
            i = ord[j]
            if !h.fine[i]
                isFine = false
                break
            end
            eLayer += h.em[i]
            c = cos_opening(h.nx[i], h.ny[i], h.nz[i], dir)
            isnan(c) && return (profile, 0.0f0, 0, false)
            x0Layer += h.x0[i] / max(abs(c), cfg.minCosAngle)
        end
        k = kEnd + 1
        isFine || break

        eCal += eLayer
        x0Layer /= Float32(nHits)
        nX0Last = x0Layer
        nX0 += x0Layer
        iLayer == inner && (nX0 *= Float32(inner + 1 - cfg.pseudoLayerAtIp))

        endPos = nX0 / bw
        endBin = min(unsafe_trunc(Int, endPos), cfg.nBins - 1)
        deltaPos = x0Layer / bw
        startPos = endPos - deltaPos
        startBin = unsafe_trunc(Int, startPos)
        for iBin in startBin:endBin
            d = 1.0f0
            if iBin == startBin
                d -= startPos - Float32(startBin)
            elseif iBin == endBin
                d -= 1.0f0 - endPos + Float32(endBin)
            end
            profile[iBin+1] += eLayer * (d / deltaPos)
        end
    end

    endBin = min(unsafe_trunc(Int, nX0 / bw), cfg.nBins)
    (endBin == 0 || eCal < eps(Float32)) && return (profile, eCal, 0, false)
    return (profile, eCal, endBin, true)
end

"""Expected profile exactly as the plugin computes it (Float32 a, float pow, double exp)."""
function expected_profile_exact(eClus::Float32, ec::Float32, cfg::PluginConfig)
    a = Float64(cfg.parameter0 + cfg.parameter1 * log(eClus / ec))
    lgA = logabsgamma(a)[1]                       # std::exp(lgamma(a)) = |Γ(a)|
    gammaA = exp(lgA)
    out = Vector{Float32}(undef, cfg.nBins)
    t = 0.0f0
    for i in 1:cfg.nBins
        t += cfg.binWidth
        v = Float64(eClus) / 2.0 * Float64((t / 2.0f0)^Float32(a - 1.0)) * exp(-Float64(t) / 2.0) *
            Float64(cfg.binWidth) / gammaA
        out[i] = Float32(v)
    end
    return out
end

"""Offset scan with early break, as in the plugin. Returns (profileStart, profileDiscrepancy, δ, ŝ)."""
function offset_scan(profile::AbstractVector{Float32}, expected::AbstractVector{Float32}, eCal::Float32,
                     endBin::Int, cfg::PluginConfig)
    minD, best, d0 = floatmax(Float32), 0, -1.0f0
    for off in 0:endBin-1
        D = 0.0f0
        for i in 0:endBin-1
            D += i < off ? profile[i+1] : abs(expected[i-off+1] - profile[i+1])
        end
        off == 0 && (d0 = D)
        if D < minD
            minD, best = D, off
        end
        D - minD > cfg.maxDifference && break
    end
    return (Float32(best) * cfg.binWidth, minD / eCal, d0 / eCal, best)
end

# =====================================================================================================================
# Fit objective
# =====================================================================================================================

"""Inputs to δ for one shower: O_i for i < N_end, E_clus (sets a and P), E_ECAL (denominator)."""
struct Shower
    O::Vector{Float32}
    eClus::Float64
    eCal::Float64
end

"""
δ(u) = (1/E_ECAL) Σ_{i<N_end} |O_i − P_i(u)|, with a = α + p1 ln E_clus. For a pure E_c fit
α = p0 − p1·u. Float64 throughout (smooth in u); same formula as the plugin.
"""
@inline function delta_alpha(sh::Shower, alpha::Float64, p1::Float64, bw::Float64)
    a = alpha + p1 * log(sh.eClus)
    lg = logabsgamma(a)[1]
    norm = sh.eClus * 0.5 * bw
    s = 0.0
    @inbounds for i in eachindex(sh.O)
        t = i * bw
        P = norm * exp((a - 1.0) * log(0.5 * t) - 0.5 * t - lg)
        s += abs(sh.O[i] - P)
    end
    return s / sh.eCal
end

@inline delta_u(sh::Shower, u::Float64, p0::Float64, p1::Float64, bw::Float64) = delta_alpha(sh, p0 - p1 * u, p1, bw)

mean_delta(shs::Vector{Shower}, u, p0, p1, bw) = mean(sh -> delta_u(sh, u, p0, p1, bw), shs)

"""Golden-section minimisation of f on [a, b]."""
function golden(f, a::Float64, b::Float64; tol = 1.0e-7, maxiter = 200)
    g = (sqrt(5.0) - 1.0) / 2.0
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ in 1:maxiter
        abs(b - a) < tol && break
        if fc < fd
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
        end
    end
    x = (a + b) / 2
    return x, f(x)
end

"""Minimum of a sampled curve: argmin plus parabola through its neighbours. Second value flags an edge minimum."""
function grid_minimum(grid::AbstractVector{Float64}, vals::AbstractVector{<:Real})
    g = argmin(vals)
    (g == 1 || g == length(grid)) && return grid[g], true
    y1, y2, y3 = Float64(vals[g-1]), Float64(vals[g]), Float64(vals[g+1])
    den = y1 - 2y2 + y3
    shift = den > 0 ? 0.5 * (y1 - y3) / den : 0.0
    return grid[g] + shift * (grid[g+1] - grid[g]), false
end

"""Matrix M[k, g] = δ_k(u_g) for bootstrap evaluation on a fixed grid."""
function delta_matrix(shs::Vector{Shower}, grid::Vector{Float64}, p0, p1, bw)
    M = Matrix{Float32}(undef, length(shs), length(grid))
    Threads.@threads for k in eachindex(shs)
        @inbounds for (g, u) in enumerate(grid)
            M[k, g] = Float32(delta_u(shs[k], u, p0, p1, bw))
        end
    end
    return M
end

"""Bootstrap replicates of û from a δ matrix: resample rows with replacement."""
function bootstrap_u(M::Matrix{Float32}, grid::Vector{Float64}, nboot::Int, rng)
    N = size(M, 1)
    us, nEdge = Vector{Float64}(undef, nboot), 0
    counts = zeros(Float32, N)
    for b in 1:nboot
        fill!(counts, 0.0f0)
        for _ in 1:N
            counts[rand(rng, 1:N)] += 1.0f0
        end
        L = vec(counts' * M) ./ N
        us[b], edge = grid_minimum(grid, L)
        nEdge += edge
    end
    return us, nEdge
end

"""Result of a one-parameter fit on a set of showers."""
struct UFit
    name::String
    n::Int
    uhat::Float64
    sigma::Float64
    edge::Bool              # minimum at an edge of the coarse scan or bootstrap grid
    coarseGrid::Vector{Float64}
    coarseL::Vector{Float64}
    boot::Vector{Float64}
end

ec_of(u) = exp(u)

function fit_u(name::String, shs::Vector{Shower}, pc::PluginConfig, fc::Dict, rng; nboot = fc["n_bootstrap"])
    p0, p1, bw = Float64(pc.parameter0), Float64(pc.parameter1), Float64(pc.binWidth)
    grid = collect(range(log(fc["ec_min_gev"]), log(fc["ec_max_gev"]); length = fc["n_coarse"]))
    L = zeros(length(grid))
    Threads.@threads for g in eachindex(grid)
        L[g] = mean_delta(shs, grid[g], p0, p1, bw)
    end
    g = argmin(L)
    edge = g == 1 || g == length(grid)
    lo, hi = grid[max(g - 1, 1)], grid[min(g + 1, length(grid))]
    uhat, _ = golden(u -> mean_delta(shs, u, p0, p1, bw), lo, hi)

    boot = Float64[]
    if nboot > 0
        hw, st = fc["fine_halfwidth"], fc["fine_step"]
        fine = collect(uhat-hw:st:uhat+hw)
        M = delta_matrix(shs, fine, p0, p1, bw)
        boot, nEdge = bootstrap_u(M, fine, nboot, rng)
        if nEdge > 0
            @warn "$name: $nEdge/$nboot bootstrap minima at the grid edge; increase fit.fine_halfwidth"
            edge = true
        end
    end
    sigma = length(boot) > 1 ? std(boot) : NaN
    return UFit(name, length(shs), uhat, sigma, edge, grid, L, boot)
end

"""F4: χ² fit of the mean normalised profile. Bins used where ≥ occ·N showers extend."""
function fit_mean_profile(shs::Vector{Shower}, pc::PluginConfig, fc::Dict, lo::Float64, hi::Float64)
    p0, p1, bw = Float64(pc.parameter0), Float64(pc.parameter1), Float64(pc.binWidth)
    nb = maximum(length(sh.O) for sh in shs)
    n = zeros(Int, nb); s = zeros(nb); s2 = zeros(nb)
    for sh in shs, i in eachindex(sh.O)
        o = sh.O[i] / sh.eCal
        n[i] += 1; s[i] += o; s2[i] += o^2
    end
    # Bins every shower leaves empty (e.g. before the extrapolated ECAL face) have zero variance and would
    # dominate χ²; F4 only uses well-populated bins carrying at least 1e-3 of the peak mean.
    occupied = n .>= fc["f4_min_bin_occupancy"] * length(shs)
    meanAll = [n[i] > 0 ? s[i] / n[i] : 0.0 for i in 1:nb]
    use = findall(occupied .& (meanAll .>= 1.0e-3 * maximum(meanAll)))
    isempty(use) && return NaN
    mo = s[use] ./ n[use]
    var = max.(s2[use] ./ n[use] .- mo .^ 2, eps()) ./ n[use]
    function chi2(u)
        mp = zeros(length(use)); cnt = zeros(Int, length(use))
        alpha = p0 - p1 * u
        for sh in shs
            a = alpha + p1 * log(sh.eClus); lg = logabsgamma(a)[1]
            for (j, i) in enumerate(use)
                i <= length(sh.O) || continue
                t = i * bw
                mp[j] += sh.eClus * 0.5 * bw * exp((a - 1) * log(0.5t) - 0.5t - lg) / sh.eCal
                cnt[j] += 1
            end
        end
        return sum((mo .- mp ./ cnt) .^ 2 ./ var)
    end
    u, _ = golden(chi2, lo, hi; tol = 1.0e-6)
    return u
end

"""Nelder–Mead in 2D (enough for the joint (α, p1) fit)."""
function nelder_mead(f, x0::Vector{Float64}, step::Vector{Float64}; tol = 1.0e-8, maxiter = 2000)
    n = length(x0)
    pts = [copy(x0)]
    for i in 1:n
        p = copy(x0); p[i] += step[i]; push!(pts, p)
    end
    vals = f.(pts)
    for _ in 1:maxiter
        o = sortperm(vals); pts, vals = pts[o], vals[o]
        abs(vals[end] - vals[1]) < tol * (abs(vals[1]) + tol) && break
        c = sum(pts[1:n]) / n
        xr = c + (c - pts[end]); fr = f(xr)
        if fr < vals[1]
            xe = c + 2(c - pts[end]); fe = f(xe)
            pts[end], vals[end] = fe < fr ? (xe, fe) : (xr, fr)
        elseif fr < vals[n]
            pts[end], vals[end] = xr, fr
        else
            xc = c + 0.5(pts[end] - c); fcv = f(xc)
            if fcv < vals[end]
                pts[end], vals[end] = xc, fcv
            else
                for i in 2:n+1
                    pts[i] = pts[1] + 0.5(pts[i] - pts[1]); vals[i] = f(pts[i])
                end
            end
        end
    end
    i = argmin(vals)
    return pts[i], vals[i]
end

"""Joint fit of (α, p1) with a = α + p1 ln E_clus, objective Σ_samples <δ> (equal weight per sample)."""
function fit_joint(sets::Vector{Vector{Shower}}, pc::PluginConfig, ustart::Float64)
    bw = Float64(pc.binWidth)
    obj(x) = sum(mean(sh -> delta_alpha(sh, x[1], x[2], bw), shs) for shs in sets)
    p0, p1 = Float64(pc.parameter0), Float64(pc.parameter1)
    x, _ = nelder_mead(obj, [p0 - p1 * ustart, p1], [0.1, 0.05])
    return x
end

# =====================================================================================================================
# Reading
# =====================================================================================================================

col(t, name) = collect(getproperty(t, Symbol(name)))

function read_config(f)
    t = LazyTree(f, "Config", ["binWidth", "nBins", "minCosAngle", "criticalEnergy", "parameter0",
                                "parameter1", "maxDifference", "pseudoLayerAtIp"])
    r = t[1]
    return PluginConfig(Float32(r.binWidth), Int(r.nBins), Float32(r.minCosAngle), Float32(r.criticalEnergy),
                        Float32(r.parameter0), Float32(r.parameter1), Float32(r.maxDifference), Int(r.pseudoLayerAtIp))
end

function read_hits(f)
    names = ["eventNumber", "clusterIndex", "isIsolated", "pseudoLayer", "hitType", "isFine", "emEnergy",
             "nCellX0", "normalX", "normalY", "normalZ", "gunWeightFraction"]
    t = LazyTree(f, "Hits", names)
    c = Dict(n => col(t, n) for n in names)
    return HitCols(UInt32.(c["eventNumber"]), Int32.(c["clusterIndex"]), Bool.(c["isIsolated"]),
                   UInt32.(c["pseudoLayer"]), Int32.(c["hitType"]), Bool.(c["isFine"]), Float32.(c["emEnergy"]),
                   Float32.(c["nCellX0"]), Float32.(c["normalX"]), Float32.(c["normalY"]), Float32.(c["normalZ"]),
                   Float32.(c["gunWeightFraction"]))
end

function read_table(f, tree, names)
    t = LazyTree(f, tree, names)
    return Dict(n => col(t, n) for n in names)
end

"""Row ranges per event number; requires (and checks) that rows of an event are contiguous."""
function event_ranges(ev::AbstractVector)
    r = Dict{UInt32,UnitRange{Int}}()
    isempty(ev) && return r
    s = 1
    for i in 2:length(ev)+1
        if i > length(ev) || ev[i] != ev[s]
            haskey(r, ev[s]) && error("rows of event $(ev[s]) are not contiguous")
            r[UInt32(ev[s])] = s:i-1
            s = i
        end
    end
    return r
end

# =====================================================================================================================
# Per-event processing
# =====================================================================================================================

mutable struct EventRec
    label::String
    file::Int
    event::UInt32
    eTrue::Float64
    cosTheta::Float64
    rEnd::Float64
    hasEndpoint::Bool
    nRoot::Int
    pdg::Int
    nClus::Int
    fFrag::Float64
    eReco::Float64
    merged::Union{Nothing,Shower}
    main::Union{Nothing,Shower}
    mainPD::Float64
    mainStart::Float64
    mainDelta::Float64
    passTruth::Bool
    passFid::Bool
    passContain::Bool
    fragClass::Char
end

mutable struct ClosureStats
    n::Int
    nPass::Int
    worstPD::Float64
    worstDelta::Float64
    worstBin::Float64
    worstEClus::Float64
    startMismatch::Int
    endBinMismatch::Int
    rows::Vector{String}
end
ClosureStats() = ClosureStats(0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, String[])

relerr(a, b, atol) = abs(a - b) / max(abs(b), atol)

function process_file!(events::Vector{EventRec}, cl::ClosureStats, path::String, ifile::Int, sample::Dict,
                       cfg::Dict, pcRef::Base.RefValue{Union{Nothing,PluginConfig}})
    label = sample["label"]
    f = ROOTFile(path)
    pc = read_config(f)
    if pcRef[] === nothing
        pcRef[] = pc
    elseif pcRef[] != pc
        error("$path: plugin Config differs from the first file of this run ($(pcRef[]) vs $pc)")
    end

    sel, clo, fc = cfg["selection"], cfg["closure"], cfg["fit"]
    wmin = Float32(sel["gun_weight_min"])
    pdgWanted = get(sample, "gun_pdg", sel["gun_pdg"])     # per-sample override, e.g. 11 for the F5 electron sample

    truth = read_table(f, "Truth", ["eventNumber", "nRootParticles", "pdg", "energy", "px", "py", "pz",
                                    "endX", "endY", "endZ"])
    prof = read_table(f, "ShowerProfiles", ["eventNumber", "clusterIndex", "nProfileBins", "measuredProfile",
                                            "profileStart", "profileDiscrepancy", "delta", "profileEndBin",
                                            "profileClusterEnergy", "eCalEnergy", "dirX", "dirY", "dirZ"])
    h = read_hits(f)

    # -1 = hit without MC information. A few untraceable hits are normal (5 of 6.5e7 in the S1+S5 baseline); a
    # sizeable fraction means MC relations are not reaching Pandora.
    nNoMC = count(==(-1.0f0), h.w)
    nNoMC > 1.0e-4 * length(h.w) && error("$path: $nNoMC/$(length(h.w)) hits with gunWeightFraction = -1 — MC relations are not reaching Pandora")
    nNoMC > 0 && @info "$path: $nNoMC hit(s) without MC information, not counted as photon hits"
    # -2 = undefined: the photon converted inside Pandora's MC PFO selection radius (r < 500 mm), so hit weights point
    # to its e+e-. Those events fail the truth selection (nRootParticles = 2); their hits are not counted as photon hits.
    h.w .= max.(h.w, 0.0f0)

    hitR = event_ranges(h.event)
    profR = event_ranges(prof["eventNumber"])
    maxClosure = clo["max_clusters_per_sample"]

    for it in eachindex(truth["eventNumber"])
        ev = UInt32(truth["eventNumber"][it])
        px, py, pz = Float64(truth["px"][it]), Float64(truth["py"][it]), Float64(truth["pz"][it])
        pmag = sqrt(px^2 + py^2 + pz^2)
        ex, ey, ez = Float64(truth["endX"][it]), Float64(truth["endY"][it]), Float64(truth["endZ"][it])
        rec = EventRec(label, ifile, ev, Float64(truth["energy"][it]), pmag > 0 ? pz / pmag : NaN, hypot(ex, ey),
                       !(ex == 0 && ey == 0 && ez == 0), Int(truth["nRootParticles"][it]), Int(truth["pdg"][it]),
                       0, NaN, 0.0, nothing, nothing, NaN, NaN, NaN, false, false, false, '-')

        hr = get(hitR, ev, 1:0)
        pr = get(profR, ev, 1:0)

        # --- closure: rebuild every dumped cluster from its own ordered hits, production direction and E_clus ---
        if clo["enabled"] && (maxClosure == 0 || cl.n < maxClosure)
            for ip in pr
                endProd = Int(prof["profileEndBin"][ip])
                (prof["nProfileBins"][ip] > 0 && endProd > 0) || continue
                ci = prof["clusterIndex"][ip]
                idx = [i for i in hr if h.cluster[i] == ci && !h.isolated[i]]
                dir = (Float32(prof["dirX"][ip]), Float32(prof["dirY"][ip]), Float32(prof["dirZ"][ip]))
                O, eCal, endBin, ok = build_profile(h, idx, dir, pc)
                eClusProd = Float32(prof["profileClusterEnergy"][ip])
                eClusHits = isempty(idx) ? 0.0f0 : sum(h.em[idx])
                cl.n += 1
                if !ok
                    push!(cl.rows, @sprintf("%s,%d,%d,%d,builder_invalid", label, ifile, ev, ci))
                    continue
                end
                P = expected_profile_exact(eClusProd, pc.criticalEnergy, pc)
                start, pd, dl, _ = offset_scan(O, P, eCal, endBin, pc)
                measured = Float32.(prof["measuredProfile"][ip])
                atol = clo["abs_tol_gev"]
                # Per-bin error relative to E_ECAL (the scale δ is normalised by). Relative to the bin's own content,
                # bins a layer only grazes amplify float-summation-order differences (~1e-7) to ~1e-4, because the
                # dump does not preserve Pandora's hit order within a layer.
                binErr = maximum(i -> abs(O[i] - measured[i]) / max(abs(eCal), Float32(atol)),
                                 1:min(endBin, length(measured)); init = 0.0)
                ePD = relerr(pd, prof["profileDiscrepancy"][ip], atol)
                eDl = relerr(dl, prof["delta"][ip], atol)
                eEc = relerr(eClusHits, eClusProd, atol)
                okStart = start == Float32(prof["profileStart"][ip])
                okEnd = endBin == endProd
                tol = clo["rel_tol"]
                pass = okStart && okEnd && ePD < tol && eDl < tol && binErr < 10tol && eEc < 10tol
                cl.nPass += pass
                cl.startMismatch += !okStart
                cl.endBinMismatch += !okEnd
                cl.worstPD = max(cl.worstPD, ePD); cl.worstDelta = max(cl.worstDelta, eDl)
                cl.worstBin = max(cl.worstBin, binErr); cl.worstEClus = max(cl.worstEClus, eEc)
                pass || push!(cl.rows, @sprintf("%s,%d,%d,%d,start %g/%g end %d/%d PD %.3g delta %.3g bin %.3g eclus %.3g",
                                                label, ifile, ev, ci, start, prof["profileStart"][ip], endBin, endProd,
                                                ePD, eDl, binErr, eEc))
            end
        end

        # --- photon energy per cluster (all hits incl. isolated, weighted by gun fraction) ---
        eByCluster = Dict{Int32,Float64}()
        for i in hr
            h.cluster[i] >= 0 || continue
            eByCluster[h.cluster[i]] = get(eByCluster, h.cluster[i], 0.0) + h.em[i] * h.w[i]
        end
        filter!(kv -> kv[2] > 0, eByCluster)
        if !isempty(eByCluster)
            mainCi, eMain = first(sort(collect(eByCluster); by = last, rev = true))
            rec.eReco = sum(values(eByCluster))
            rec.nClus = length(eByCluster)
            rec.fFrag = 1 - eMain / rec.eReco

            # main cluster, production inputs
            ipMain = findfirst(ip -> prof["clusterIndex"][ip] == mainCi, pr)
            mainDir = nothing
            if ipMain !== nothing
                ip = pr[ipMain]
                mainDir = (Float32(prof["dirX"][ip]), Float32(prof["dirY"][ip]), Float32(prof["dirZ"][ip]))
                endProd = Int(prof["profileEndBin"][ip])
                if prof["nProfileBins"][ip] > 0 && endProd > 0 && prof["eCalEnergy"][ip] > 0
                    m = Float32.(prof["measuredProfile"][ip])
                    rec.main = Shower(m[1:endProd], Float64(prof["profileClusterEnergy"][ip]), Float64(prof["eCalEnergy"][ip]))
                    rec.mainPD = prof["profileDiscrepancy"][ip]
                    rec.mainStart = prof["profileStart"][ip]
                    rec.mainDelta = prof["delta"][ip]
                end
            end

            # merged hit set: ordered (non-isolated) photon hits in any cluster
            idx = [i for i in hr if h.cluster[i] >= 0 && !h.isolated[i] && h.w[i] >= wmin]
            dir = if fc["direction"] == "truth"
                pmag > 0 ? (Float32(px / pmag), Float32(py / pmag), Float32(pz / pmag)) : nothing
            elseif fc["direction"] == "mainfit"
                mainDir
            else
                error("fit.direction must be \"truth\" or \"mainfit\"")
            end
            if dir !== nothing && !isempty(idx)
                O, eCal, endBin, ok = build_profile(h, idx, dir, pc)
                eClus = sum(Float64, h.em[idx])
                ok && (rec.merged = Shower(O[1:endBin], eClus, Float64(eCal)))
            end
        end

        rec.passTruth = rec.pdg == pdgWanted &&
                        (!sel["require_single_root"] || rec.nRoot == 1) &&
                        (!sel["require_endpoint"] || rec.hasEndpoint) &&
                        (!rec.hasEndpoint || rec.rEnd >= sel["min_conversion_radius_mm"])
        rec.passFid = abs(rec.cosTheta) <= sel["max_abs_costheta"]
        edges = fc["fragmentation_edges"]
        rec.fragClass = isnan(rec.fFrag) ? '-' : rec.fFrag < edges[1] ? 'A' : rec.fFrag < edges[2] ? 'B' : 'C'
        push!(events, rec)
    end
    return nothing
end

"""Containment cut per sample: keep E_reco/E_true at or above the given percentile (of truth+fiducial events)."""
function apply_containment!(evs::Vector{EventRec}, pct::Float64)
    cand = [e for e in evs if e.passTruth && e.passFid && e.merged !== nothing]
    isempty(cand) && return NaN
    thr = pct > 0 ? quantile([e.eReco / e.eTrue for e in cand], pct / 100) : -Inf
    for e in evs
        e.passContain = e.merged !== nothing && e.eReco / e.eTrue >= thr
    end
    return thr
end

selected(e::EventRec) = e.passTruth && e.passFid && e.passContain && e.merged !== nothing

# =====================================================================================================================
# Output helpers
# =====================================================================================================================

function write_csv(path, header::Vector{String}, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for r in rows
            println(io, join(r, ","))
        end
    end
end

fmt(x::AbstractFloat) = isnan(x) ? "nan" : @sprintf("%.8g", x)
fmt(x) = string(x)

chi2_pvalue(chi2, ndf) = ndf > 0 ? gamma_inc(ndf / 2, chi2 / 2)[2] : NaN

# =====================================================================================================================
# Main
# =====================================================================================================================

function main(cfgPath::String)
    cfg = merge_defaults(TOML.parsefile(cfgPath), DEFAULTS)
    sel, fc, clo = cfg["selection"], cfg["fit"], cfg["closure"]
    outdir = cfg["output_dir"]
    mkpath(outdir)
    rng = MersenneTwister(fc["seed"])
    log_io = open(joinpath(outdir, "summary.txt"), "w")
    say(args...) = (println(args...); println(log_io, args...); flush(log_io))

    samples = cfg["sample"]
    isempty(samples) && error("no [[sample]] entries in $cfgPath")
    say("ecFit.jl — ", length(samples), " samples, direction = ", fc["direction"], ", threads = ", Threads.nthreads())

    pcRef = Ref{Union{Nothing,PluginConfig}}(nothing)
    allEvents = Dict{String,Vector{EventRec}}()
    closure = ClosureStats()

    # ---------------- read ----------------
    for s in samples
        label = s["label"]
        files = reduce(vcat, expand_glob.(s["files"]))
        isempty(files) && error("sample $label: no files match $(s["files"])")
        evs = EventRec[]
        for (ifile, path) in enumerate(files)
            process_file!(evs, closure, path, ifile, s, cfg, pcRef)
        end
        thr = apply_containment!(evs, Float64(sel["containment_percentile"]))
        allEvents[label] = evs
        nT = count(e -> e.passTruth, evs); nF = count(e -> e.passTruth && e.passFid, evs)
        say(@sprintf("  %-14s files %4d  events %7d  truth %7d  fiducial %7d  selected %7d  (E_reco/E_true cut %.4f)",
                     label, length(files), length(evs), nT, nF, count(selected, evs), thr))
    end
    pc = pcRef[]
    say(@sprintf("Plugin config: binWidth %.3g X0, nBins %d, minCos %.3g, Ec %.6g GeV, p0 %.4g, p1 %.4g, maxDiff %.3g, layerAtIp %d",
                 pc.binWidth, pc.nBins, pc.minCosAngle, pc.criticalEnergy, pc.parameter0, pc.parameter1,
                 pc.maxDifference, pc.pseudoLayerAtIp))

    # ---------------- closure ----------------
    if clo["enabled"]
        frac = closure.n > 0 ? closure.nPass / closure.n : NaN
        say(@sprintf("Closure: %d/%d clusters pass (%.4f%%); start mismatches %d, N_end mismatches %d; worst rel. err: PD %.2g, δ %.2g, bin %.2g, E_clus(hits) %.2g",
                     closure.nPass, closure.n, 100frac, closure.startMismatch, closure.endBinMismatch,
                     closure.worstPD, closure.worstDelta, closure.worstBin, closure.worstEClus))
        write_csv(joinpath(outdir, "closure_failures.csv"), ["label", "file", "event", "cluster", "detail"],
                  [split(r, ","; limit = 5) for r in closure.rows])
        if !(frac >= clo["require_pass_fraction"])
            say("ABORT: closure pass fraction below closure.require_pass_fraction — see closure_failures.csv")
            close(log_io)
            exit(2)
        end
    end

    ec0 = Float64(pc.criticalEnergy)
    p0, p1, bw = Float64(pc.parameter0), Float64(pc.parameter1), Float64(pc.binWidth)

    # ---------------- F1 + F4 per sample ----------------
    results = Vector{Vector{String}}()
    scans = Vector{Vector{String}}()
    boots = Vector{Vector{String}}()
    sampleFits = Dict{String,UFit}()
    sampleShowers = Dict{String,Vector{Shower}}()
    say("\nF1 (merged showers, δ at s = 0) and F4 (mean-profile χ²):")
    for s in samples
        label = s["label"]
        shs = [e.merged for e in allEvents[label] if selected(e)]
        if length(shs) < 10
            say("  $label: only $(length(shs)) selected showers — skipped")
            continue
        end
        fit = fit_u(label, shs, pc, fc, rng)
        uF4 = fit_mean_profile(shs, pc, fc, fit.uhat - 1.5, fit.uhat + 1.5)
        sampleFits[label] = fit; sampleShowers[label] = shs
        d0 = mean_delta(shs, log(ec0), p0, p1, bw); d1 = mean_delta(shs, fit.uhat, p0, p1, bw)
        say(@sprintf("  %-14s N %6d  Ec = %.5g GeV  (u = %.4f ± %.4f)%s   <δ> default %.4f → fitted %.4f   F4: Ec = %.5g (Δu = %+.4f)",
                     label, fit.n, ec_of(fit.uhat), fit.uhat, fit.sigma, fit.edge ? " [EDGE]" : "", d0, d1,
                     ec_of(uF4), uF4 - fit.uhat))
        push!(results, fmt.(["F1", label, s["group"], s["energy_gev"], fit.n, fit.uhat, fit.sigma, ec_of(fit.uhat),
                             ec_of(fit.uhat - fit.sigma), ec_of(fit.uhat + fit.sigma), fit.edge, uF4, d0, d1]))
        for (u, L) in zip(fit.coarseGrid, fit.coarseL)
            push!(scans, fmt.(["sample", label, u, ec_of(u), L]))
        end
        for (b, u) in enumerate(fit.boot)
            push!(boots, fmt.(["sample", label, b, u]))
        end
    end

    # ---------------- F2 per group ----------------
    groups = unique(s["group"] for s in samples)
    groupU = Dict{String,Tuple{Float64,Float64}}()
    adopted = nothing
    say("\nF2 (per group):")
    for grp in groups
        labs = [s["label"] for s in samples if s["group"] == grp && haskey(sampleFits, s["label"])]
        isempty(labs) && continue
        us = [sampleFits[l].uhat for l in labs]; σs = [sampleFits[l].sigma for l in labs]
        w = 1 ./ σs .^ 2
        uw = sum(w .* us) / sum(w); σw = 1 / sqrt(sum(w))
        chi2 = sum(w .* (us .- uw) .^ 2); ndf = length(us) - 1; pval = chi2_pvalue(chi2, ndf)
        sets = [sampleShowers[l] for l in labs]

        # combined, equal weight per sample: Σ_j <δ>_j
        objC(u) = sum(mean_delta(shs, u, p0, p1, bw) for shs in sets)
        uc, _ = golden(objC, minimum(us) - 0.5, maximum(us) + 0.5)
        hw, st = fc["fine_halfwidth"], fc["fine_step"]
        fine = collect(minimum(us)-hw:st:maximum(us)+hw)
        if length(fine) > 1500
            fine = collect(range(first(fine), last(fine); length = 1500))
        end
        Ms = [delta_matrix(shs, fine, p0, p1, bw) for shs in sets]
        bootC = Float64[]
        for _ in 1:fc["n_bootstrap"]
            L = zeros(length(fine))
            for M in Ms
                N = size(M, 1); c = zeros(Float32, N)
                for _ in 1:N
                    c[rand(rng, 1:N)] += 1.0f0
                end
                L .+= vec(c' * M) ./ N
            end
            push!(bootC, grid_minimum(fine, L)[1])
        end
        σc = std(bootC)
        groupU[grp] = (uc, σc)
        for (b, u) in enumerate(bootC)
            push!(boots, fmt.(["group", grp, b, u]))
        end
        for u in fine[1:10:end]
            push!(scans, fmt.(["group", grp, u, ec_of(u), objC(u)]))
        end

        # joint (α, p1): needs at least two distinct energies to separate α from p1
        energies = unique(Float64(s["energy_gev"]) for s in samples if s["label"] in labs)
        x, σα, σp1, ecJoint = [NaN, NaN], NaN, NaN, NaN
        if length(energies) >= 2
            x = fit_joint(sets, pc, uc)
            jb = Vector{Vector{Float64}}()
            for _ in 1:fc["n_bootstrap_joint"]
                push!(jb, fit_joint([shs[rand(rng, 1:length(shs), length(shs))] for shs in sets], pc, uc))
            end
            σα = length(jb) > 1 ? std(first.(jb)) : NaN
            σp1 = length(jb) > 1 ? std(last.(jb)) : NaN
            ecJoint = exp((p0 - x[1]) / x[2])
        end

        spread = maximum(us) - minimum(us)
        constant = length(us) < 2 || pval > fc["constancy_pvalue"] || spread < fc["constancy_min_spread_u"]
        say(@sprintf("  group %-8s %d samples: χ²/ndf = %.2f/%d (p = %.3g), spread in u = %.4f → %s", grp, length(us),
                     chi2, ndf, pval, spread, constant ? "consistent with constant E_c" : "TREND with energy"))
        say(@sprintf("    combined (equal weight):   Ec = %.5g GeV  (u = %.4f ± %.4f)", ec_of(uc), uc, σc))
        say(@sprintf("    inverse-variance mean:     Ec = %.5g GeV  (u = %.4f ± %.4f)  [systematic variant]", ec_of(uw), uw, σw))
        if isnan(ecJoint)
            say("    joint: skipped (needs ≥ 2 energies)")
        else
            say(@sprintf("    joint: α = %.4f ± %.4f, p1 = %.4f ± %.4f  → effective Ec (p0 = %.3g) = %.5g GeV",
                         x[1], σα, x[2], σp1, p0, ecJoint))
        end
        push!(results, fmt.(["F2_combined", grp, grp, NaN, sum(length, sets), uc, σc, ec_of(uc), ec_of(uc - σc),
                             ec_of(uc + σc), false, NaN, NaN, NaN]))
        push!(results, fmt.(["F2_invvar", grp, grp, NaN, sum(length, sets), uw, σw, ec_of(uw), ec_of(uw - σw),
                             ec_of(uw + σw), false, chi2, ndf, pval]))
        push!(results, fmt.(["F2_joint", grp, grp, NaN, sum(length, sets), x[1], σα, ecJoint, x[2], σp1, false,
                             NaN, NaN, NaN]))
        if grp == cfg["primary_group"]
            adopted = constant ? (:ec, ec_of(uc), σc) : (:joint, ecJoint, x[2])
        end
    end

    if haskey(groupU, cfg["primary_group"])
        uP, σP = groupU[cfg["primary_group"]]
        for (grp, (u, σ)) in groupU
            grp == cfg["primary_group"] && continue
            say(@sprintf("  Ec(%s)/Ec(%s) = %.4f ± %.4f   (electron closure F5 expects e ≈ 2.718)", grp,
                         cfg["primary_group"], exp(u - uP), exp(u - uP) * sqrt(σ^2 + σP^2)))
        end
    end

    # ---------------- F3: main cluster by fragmentation class ----------------
    say("\nF3 (main cluster, production inputs, by fragmentation class):")
    for s in samples
        label = s["label"]
        haskey(sampleFits, label) || continue
        sel_ev = [e for e in allEvents[label] if selected(e) && e.main !== nothing]
        line = @sprintf("  %-14s", label)
        for cls in ('A', 'B', 'C')
            shs = [e.main for e in sel_ev if e.fragClass == cls]
            if length(shs) < fc["min_events_per_class"]
                line *= @sprintf("  %c: n=%-6d —                ", cls, length(shs))
                continue
            end
            fit = fit_u("$label/$cls", shs, pc, fc, rng; nboot = max(fc["n_bootstrap"] ÷ 5, 20))
            line *= @sprintf("  %c: n=%-6d Δu = %+.4f ± %.4f", cls, fit.n, fit.uhat - sampleFits[label].uhat, fit.sigma)
            push!(results, fmt.(["F3_" * cls, label, s["group"], s["energy_gev"], fit.n, fit.uhat, fit.sigma,
                                 ec_of(fit.uhat), ec_of(fit.uhat - fit.sigma), ec_of(fit.uhat + fit.sigma), fit.edge,
                                 NaN, NaN, NaN]))
        end
        say(line)
    end

    # ---------------- mean profiles (for overlays) ----------------
    mp = Vector{Vector{String}}()
    # Adopted parameters as a = α + p1·ln E_clus. Constant-E_c adoption: α = p0 − p1·ln E_c with p1 unchanged.
    # Joint adoption: (α, p1) from the joint fit, recovered from the effective E_c defined at p0.
    alphaDef = p0 - p1 * log(ec0)
    alphaFit, p1Fit = if adopted === nothing
        alphaDef, p1
    elseif adopted[1] == :ec
        p0 - p1 * log(adopted[2]), p1
    else
        p0 - adopted[3] * log(adopted[2]), adopted[3]
    end
    for (label, shs) in sampleShowers
        nb = maximum(length(sh.O) for sh in shs)
        for i in 1:nb
            o = [sh.O[i] / sh.eCal for sh in shs if i <= length(sh.O)]
            pDef = Float64[]; pFit = Float64[]
            for sh in shs
                i <= length(sh.O) || continue
                for (α, pp1, v) in ((alphaDef, p1, pDef), (alphaFit, p1Fit, pFit))
                    a = α + pp1 * log(sh.eClus); lg = logabsgamma(a)[1]; t = i * bw
                    push!(v, sh.eClus * 0.5 * bw * exp((a - 1) * log(0.5t) - 0.5t - lg) / sh.eCal)
                end
            end
            push!(mp, fmt.([label, i, i * bw, length(o), mean(o), std(o), mean(pDef), mean(pFit)]))
        end
    end

    # ---------------- per-event table ----------------
    evRows = Vector{Vector{String}}()
    for s in samples, e in allEvents[s["label"]]
        dDef = e.merged === nothing ? NaN : delta_u(e.merged, log(ec0), p0, p1, bw)
        dFit = e.merged === nothing ? NaN : delta_alpha(e.merged, alphaFit, p1Fit, bw)
        push!(evRows, fmt.([e.label, e.file, e.event, e.eTrue, e.cosTheta, e.rEnd, e.hasEndpoint, e.nRoot, e.pdg,
                            e.nClus, e.fFrag, e.eReco, e.merged !== nothing, e.main !== nothing, e.passTruth, e.passFid,
                            e.passContain, selected(e), e.fragClass, dDef, dFit, e.mainPD, e.mainStart, e.mainDelta]))
    end

    write_csv(joinpath(outdir, "results.csv"),
              ["fit", "name", "group", "energy_gev", "n", "u_or_alpha", "sigma", "ec_gev", "ec_lo_or_p1", "ec_hi_or_sigma_p1",
               "edge", "uF4_or_chi2", "meanDelta_default_or_ndf", "meanDelta_fit_or_p"], results)
    write_csv(joinpath(outdir, "scans.csv"), ["scope", "name", "u", "ec_gev", "L"], scans)
    write_csv(joinpath(outdir, "bootstrap.csv"), ["scope", "name", "replicate", "u"], boots)
    write_csv(joinpath(outdir, "mean_profiles.csv"),
              ["sample", "bin", "t_upper_x0", "n", "mean_o", "std_o", "mean_p_default", "mean_p_fit"], mp)
    write_csv(joinpath(outdir, "events.csv"),
              ["label", "file", "event", "e_true", "cos_theta", "r_end_mm", "has_endpoint", "n_root", "pdg", "n_clus",
               "f_frag", "e_reco", "merged_ok", "main_ok", "pass_truth", "pass_fid", "pass_contain", "selected",
               "frag_class", "delta_merged_default", "delta_merged_fit", "pd_main", "start_main", "delta_main"], evRows)

    # ---------------- decision ----------------
    say("")
    if adopted === nothing
        say("No fit for primary group \"$(cfg["primary_group"])\".")
    elseif adopted[1] == :ec
        say(@sprintf("ADOPT (primary group %s, constant E_c): Ec = %.5g GeV, σ(ln Ec) = %.4f",
                     cfg["primary_group"], adopted[2], adopted[3]))
        say(@sprintf("  <LongProfileCriticalEnergy>%.6g</LongProfileCriticalEnergy>", adopted[2]))
    else
        say(@sprintf("ADOPT (primary group %s, energy trend): joint fit, p1 = %.5g, effective Ec = %.5g GeV",
                     cfg["primary_group"], adopted[3], adopted[2]))
        say(@sprintf("  <LongProfileCriticalEnergy>%.6g</LongProfileCriticalEnergy>", adopted[2]))
        say(@sprintf("  <LongProfileParameter1>%.6g</LongProfileParameter1>", adopted[3]))
    end
    if adopted !== nothing
        say("  <δ> per sample, default → adopted parameters:")
        for s in samples
            haskey(sampleShowers, s["label"]) || continue
            shs = sampleShowers[s["label"]]
            dDef = mean(sh -> delta_alpha(sh, alphaDef, p1, bw), shs)
            dAdp = mean(sh -> delta_alpha(sh, alphaFit, p1Fit, bw), shs)
            say(@sprintf("    %-14s %.4f → %.4f  (%+.1f%%)", s["label"], dDef, dAdp, 100 * (dAdp / dDef - 1)))
        end
    end
    say("Outputs in $(abspath(outdir))")
    close(log_io)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || (println("usage: julia --project=<env> ecFit.jl ecfit.toml"); exit(1))
    main(ARGS[1])
end
