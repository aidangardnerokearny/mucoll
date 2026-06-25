 using UnROOT 
 using Plots
 using StatsBase
 using Statistics
 using SpecialFunctions
 using LaTeXStrings

 # --- Arguments ---
 file_nom = length(ARGS) >= ARGS[1] : "shower_profiles.root"
 file_vac = length(ARGS) >= ARGS[2] : "shower_profiles_vacuum.root"
 output_file = length(ARGS) >= ARGS[3] : "vacComp.pdf"
 energy_min = length(ARGS) >= ARGS[4] : 0.1
 energy_max = length(ARGS) >= ARGS[5] : 1e9


 # --- Open Files ---
 println("Opening $file_nom...")
 f_nom = ROOTFile(file_nom)
 println("Opening $file_vac...")
 f_vac = ROOTFile(file_vac)


