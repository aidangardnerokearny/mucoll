# E_c fit (S1 baseline, 16 Nov 2026)

This is the fit of `LongProfileCriticalEnergy` described in `jobs/fit_plan.md`. The inputs are the dumps from the S1 baseline
reconstruction (vacuum solenoid, single photons, E_c = 80 MeV) in `jobs/output/ec_fit/S1/*/reco_baseline/`.

```
julia -t 16 Ecfit.jl ecfit.toml            # primary: truth direction             → files in this directory
julia -t 16 Ecfit.jl ecfit_mainfit.toml    # systematic: main-cluster fit direction → direction_mainfit/
julia plot_ecfit.jl [dir]                  # PDFs from a fit's CSVs
```

The run uses the global Julia environment, which needs UnROOT, SpecialFunctions, CSV, DataFrames and CairoMakie.

| File | Content |
|---|---|
| `summary.txt` | Full printed result: selection, closure, F1–F4 and the decision |
| `results.csv` | One row per fit (F1 per energy, F2 combined / inverse-variance / joint, F3 per class) |
| `scans.csv`, `bootstrap.csv` | L(u) scans and bootstrap replicas |
| `mean_profiles.csv` | ⟨O_i⟩/E_ECAL with P_i at the default and the adopted parameters |
| `events.csv` | Per-event selection flags, fragmentation, δ at the default and the adopted parameters, and main-cluster PD/start |
| `closure_failures.csv` | Closure failures (empty: all 93 420 clusters pass) |
| `ecfit_*.pdf` | Scans, E_c vs energy, mean profiles with residuals, F3 by fragmentation class |
| `fit_run.log`, `direction_mainfit_run.log` | Console output of the two runs |

## Changes made to Ecfit.jl for this run

- `gunWeightFraction` −2 (photon converted at r < 500 mm) and −1 (no MC information; 5 of 6.5×10⁷ hits) no longer
  abort the run. Such hits are not counted as photon hits. The script still aborts if more than 10⁻⁴ of a file's hits are −1.
- The per-bin closure check is now normalised by E_ECAL rather than by the bin's own content. On bins that a layer
  only grazes, the old check amplified 10⁻⁷ float-ordering differences above tolerance. profileStart and N_end must still
  match exactly, and PD and δ to 10⁻⁵.
- Bug fix: when the joint (α, p1) fit is adopted, `mean_profiles.csv` and `events.csv` now use the adopted (α, p1).
  Before the fix they silently used the default 80 MeV.
- `summary.txt` now lists ⟨δ⟩ per energy at the default and at the adopted parameters.
