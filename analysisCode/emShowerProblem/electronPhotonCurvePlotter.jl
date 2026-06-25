using Plots
using SpecialFunctions
using LaTeXStrings

# --- Parameters ---
energies = [1.0, 5.0, 25.0, 100.0]
C_i = [-0.5, 0.5]
b = 0.5
E_c = 0.08

expected_profile = []

# --- Iterate over all energies and both values of C_i ---
for energy in energies
    for i in C_i
        profile = []
        a =  1 + i + b*log(energy/E_c)
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

for (subplot_idx, data_idx) in enumerate(1:2:length(expected_profile))
    y1 = expected_profile[data_idx]      # C_i = -0.5 (electron)
    y2 = expected_profile[data_idx + 1]  # C_i =  0.5 (photon)
    energy = energies[subplot_idx]
    plot!(p[subplot_idx], ts, y1, label = L"C_e=-0.5",      linewidth = 2, xlabel=L"X_0", ylabel=L"E \; \mathrm{ [GeV]}")
    plot!(p[subplot_idx], ts, y2, label = L"C_\gamma=+0.5", linewidth = 2)
    title!(p[subplot_idx], "E = $energy GeV")
end

savefig("showerProblemPlots/electron_photon_curve.pdf")
