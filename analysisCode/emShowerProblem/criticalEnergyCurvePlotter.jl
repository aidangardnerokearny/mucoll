using Plots
using SpecialFunctions
using LaTeXStrings

# --- Parameters ---
energies = [15.0, 20.0, 25.0, 30.0]
C_i = -0.5
b = 0.5
E_c = [0.005, 0.01, 0.08, 0.15]

expected_profile = []

# --- Iterate over all energies and all values of E_c ---
for energy in energies
    for i in E_c
        profile = []
        a =  1 + C_i + b*log(energy/i)
        gamma_a = exp(SpecialFunctions.loggamma(a))

        t = 0.0
        while t <= 49.5 
            t += 0.5
            push!(profile, energy * b * ((b*t) ^ (a-1)) * exp(-b * t) * 0.5 / gamma_a)
        end
        push!(expected_profile, profile)
    end
end

println("Data shape: ", size(expected_profile))
println("Single function size: ", length(expected_profile[1]))

# --- Plot Defaults ---
default(
        grid = false,
        framestyle = :box,
        legend_background_color = :transparent,
        legend_foreground_color = :transparent,
        size = (1000, 1000),
        legend = :topright,
        legendfontsize = 12,
)


# --- Plot Curves ---
ts = range(0, 50, length=length(expected_profile[1]))
p = plot(layout=(2,2))

for (subplot_idx, data_idx) in enumerate(1:4:length(expected_profile))
    y1 = expected_profile[data_idx]      # C_i = 0.04
    y2 = expected_profile[data_idx + 1]  # C_i = 0.06 
    y3 = expected_profile[data_idx + 2]  # E_c = 0.08
    y4 = expected_profile[data_idx + 3]  # E_c = 0.1
    energy = energies[subplot_idx]
    plot!(p[subplot_idx], ts, y1, label = L"E_c=0.005", linewidth = 2, xlabel=L"X_0", ylabel=L"E \; \mathrm{ [GeV]}")
    plot!(p[subplot_idx], ts, y2, label = L"E_c=0.01", linewidth = 2)
    plot!(p[subplot_idx], ts, y3, label = L"E_c=0.08", linewidth = 2)
    plot!(p[subplot_idx], ts, y4, label = L"E_c=0.15", linewidth = 2)
    title!(p[subplot_idx], "E = $energy GeV")
end

savefig("showerProblemPlots/critical_energy_curve.pdf")
