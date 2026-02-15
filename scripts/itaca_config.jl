#!/usr/bin/env julia

"""
ITACA Configuration & Basic Calculations
=========================================
Prints ITACA default parameters and computes basic derived values.
"""

include("itaca_f.jl")

function main()
    println("\n" * "=" ^ 60)
    println("          ITACA Detector Properties")
    println("=" ^ 60)

    # =================================================================
    # Default ITACA Configuration
    # =================================================================
    println("\n--- Default ITACA Configuration ---")
    println("\nGeometry:")
    @printf("  Fiducial length:    L_fid = %.1f m\n", L_fid_default)
    @printf("  Fiducial diameter:  D_fid = %.1f m\n", D_fid_default)
    println("\nElectric Field & Drift:")
    @printf("  Electric field:     E = %.0f V/cm\n", E_default)
    @printf("  Ion drift velocity: v_drift = %.0f cm/s\n", v_drift_default)
    println("\nVoltages:")
    @printf("  Cathode voltage:    V_c = %.1f kV\n", V_cathode)
    @printf("  EL voltage:         V_EL = %.1f kV\n", V_EL)
    println("\nGas Properties:")
    @printf("  Xenon density (15 bar): ρ = %.1f kg/m³\n", ρ_Xe_15bar)

    # =================================================================
    # Compute ITACA Default Values
    # =================================================================
    println("\n" * "-" ^ 60)
    println("--- ITACA Default Configuration Results ---")
    println("-" ^ 60)

    # Fiducial mass
    mass_kg = fiducial_mass_kg(D_fid_default, L_fid_default)
    mass_ton = fiducial_mass_ton(D_fid_default, L_fid_default)
    println("\nFiducial Mass:")
    @printf("  Mass = %.1f kg = %.3f ton\n", mass_kg, mass_ton)

    # Ion diffusion
    σ_mm = transverse_diffusion_mm(L_fid_default * 1000.0, E_default)
    println("\nIon Transverse Diffusion:")
    @printf("  σ = %.3f mm (at L = %.1f m, E = %.0f V/cm)\n", σ_mm, L_fid_default, E_default)

    # Electrode voltages
    V1 = voltage_electrode_1(V_cathode, E_default, L_fid_default)
    V2 = voltage_electrode_2(V_cathode, E_default, L_fid_default)
    println("\nFAT-GEM Electrode Voltages:")
    @printf("  V_electrode_1 = V_c + E × L = %.1f + %.0f × %.1f = %.1f kV\n",
            V_cathode, E_default, L_fid_default * 100, V1)
    @printf("  V_electrode_2 = V_electrode_1 + V_EL = %.1f + %.1f = %.1f kV\n",
            V1, V_EL, V2)

    # Drift times
    ion_drift_s = ion_drift_time_s(L_fid_default * 100.0)
    electron_drift_ms = electron_drift_time_ms(L_fid_default * 100.0)
    println("\nDrift Times:")
    @printf("  Ion drift time:      t_ion = %.1f s (full drift)\n", ion_drift_s)
    @printf("  Electron drift time: t_e = %.2f ms (full drift)\n", electron_drift_ms)

    # ASME Pressure Vessel Thickness (for default diameter)
    P_bar = 15.0
    S_MPa = 137.0  # Stainless steel 316L
    thicknesses = pressure_vessel_thicknesses(P_bar, D_fid_default; S_MPa=S_MPa, E_w=1.0)
    println("\nASME Pressure Vessel Thickness (Section VIII, Div. 1):")
    @printf("  Operating pressure:    P = %.0f bar (%.1f MPa)\n", P_bar, P_bar * 0.1)
    @printf("  Inner diameter:        D_i = %.1f m\n", D_fid_default)
    @printf("  Material:              SS 316L (S = %.0f MPa)\n", S_MPa)
    @printf("  Joint efficiency:      E_w = 1.0 (full radiography)\n")
    @printf("  Cylinder wall:         t_cyl = %.2f mm\n", thicknesses.t_cylinder_mm)
    @printf("  Hemisphere wall:       t_hem = %.2f mm\n", thicknesses.t_hemisphere_mm)

    # Cathode mesh transparency
    f = mesh_transparency(w_mesh_default, M_mesh_default)
    println("\nCathode Mesh:")
    @printf("  Wire diameter:   w = %.0f μm\n", w_mesh_default * 1e6)
    @printf("  Pitch:           M = %.0f mm\n", M_mesh_default * 1e3)
    @printf("  Transparency:    f = %.2f (%.0f%%)\n", f, f * 100)

    println("\nDone!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
