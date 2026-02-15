#!/usr/bin/env julia

"""
ITACA Plot Generation
=====================
Generates all plots for the ITACA paper.
"""

include("itaca_f.jl")

function main()
    println("\n" * "=" ^ 60)
    println("          ITACA Plot Generation")
    println("=" ^ 60)

    # Create output directory
    script_dir = @__DIR__
    output_dir = joinpath(script_dir, "itacaPlots")
    mkpath(output_dir)
    println("\nOutput directory: $output_dir\n")

    # =================================================================
    # Plot 1: Fiducial Mass vs Diameter
    # =================================================================
    println("  [1/4] Fiducial mass vs diameter...")
    p1 = plot_fiducial_mass_vs_diameter(; L_values=[1.0, 1.5, 2.0], D_min=1.0, D_max=4.0)
    savefig(p1, joinpath(output_dir, "fiducial_mass_vs_diameter.png"))
    println("        -> fiducial_mass_vs_diameter.png")

    # =================================================================
    # Plot 2: Ion Diffusion vs Drift Length (multiple E fields)
    # =================================================================
    println("  [2/4] Ion diffusion vs drift length...")
    p2 = plot_ion_diffusion_vs_length(; E_values=[200.0, 300.0, 400.0], L_max=1.5)
    savefig(p2, joinpath(output_dir, "ion_diffusion_vs_length.png"))
    println("        -> ion_diffusion_vs_length.png")

    # =================================================================
    # Plot 3: Electrode Voltage vs Drift Length
    # =================================================================
    println("  [3/4] Electrode voltage vs drift length...")
    p3 = plot_electrode_voltage_vs_length(; E_values=[200.0, 300.0, 400.0], V_c=V_cathode, L_max=1.5)
    savefig(p3, joinpath(output_dir, "electrode_voltage_vs_length.png"))
    println("        -> electrode_voltage_vs_length.png")

    # =================================================================
    # Plot 4: Electron vs Ion Diffusion
    # =================================================================
    println("  [4/4] Electron vs ion diffusion...")
    p4 = plot_electron_ion_diffusion(; E=300.0, P_bar=15.0, L_max_cm=150.0)
    savefig(p4, joinpath(output_dir, "electron_ion_diffusion.png"))
    println("        -> electron_ion_diffusion.png")

    # =================================================================
    # Summary
    # =================================================================
    println("\n" * "-" ^ 40)
    println("Generated 4 plots in: $output_dir")
    println("-" ^ 40)
    println("  1. fiducial_mass_vs_diameter.png")
    println("  2. ion_diffusion_vs_length.png")
    println("  3. electrode_voltage_vs_length.png")
    println("  4. electron_ion_diffusion.png")

    println("\nDone!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
