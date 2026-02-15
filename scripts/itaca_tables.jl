#!/usr/bin/env julia

"""
ITACA Parameter Tables
======================
Generates parameter sweep tables for design exploration.
"""

include("itaca_f.jl")

function main()
    println("\n" * "=" ^ 60)
    println("          ITACA Parameter Tables")
    println("=" ^ 60)

    # =================================================================
    # Fiducial Mass Table
    # =================================================================
    println("\n--- Fiducial Masses (ρ = $(ρ_Xe_15bar) kg/m³ at 15 bar) ---")
    println("\nD_fid (m)  |  L_fid (m)  |  Mass (kg)  |  Mass (ton)")
    println("-" ^ 55)

    for D in [1.0, 2.0, 2.6, 3.2, 4.0]
        for L in [1.0, 1.5, 2.0]
            m_kg = fiducial_mass_kg(D, L)
            m_ton = fiducial_mass_ton(D, L)
            @printf("   %.1f      |    %.1f      |  %8.1f   |   %.3f\n", D, L, m_kg, m_ton)
        end
    end

    # =================================================================
    # Ion Transverse Diffusion Table
    # =================================================================
    println("\n--- Ion Transverse Diffusion (T = 300 K) ---")
    println("\nL (m)  |  E (V/cm)  |  σ (mm)")
    println("-" ^ 40)

    for L in [0.5, 1.0, 1.5]
        for E in [200.0, 300.0, 400.0]
            σ = transverse_diffusion_mm(L * 1000.0, E)
            @printf(" %.1f   |    %3.0f     |  %.3f\n", L, E, σ)
        end
    end

    # =================================================================
    # FAT-GEM Electrode Voltages Table
    # =================================================================
    println("\n--- FAT-GEM Electrode Voltages (V_c = $(V_cathode) kV, V_EL = $(V_EL) kV) ---")
    println("\nL (m)  |  E (V/cm)  |  V_elec1 (kV)  |  V_elec2 (kV)")
    println("-" ^ 55)

    for L in [0.5, 1.0, 1.5]
        for E in [200.0, 300.0, 400.0]
            V1 = voltage_electrode_1(V_cathode, E, L)
            V2 = voltage_electrode_2(V_cathode, E, L)
            @printf(" %.1f   |    %3.0f     |     %5.1f      |     %5.1f\n", L, E, V1, V2)
        end
    end

    # =================================================================
    # ASME Pressure Vessel Wall Thicknesses Table
    # =================================================================
    println("\n--- ASME Pressure Vessel Wall Thicknesses (P = 15 bar, SS 316L) ---")
    println("\nD_i (m)  |  t_cyl (mm)  |  t_hem (mm)")
    println("-" ^ 42)

    for D in [1.0, 2.0, 2.6, 3.2, 4.0]
        th = pressure_vessel_thicknesses(15.0, D; S_MPa=137.0, E_w=1.0)
        @printf("  %.1f     |    %5.2f     |    %5.2f\n", D, th.t_cylinder_mm, th.t_hemisphere_mm)
    end

    # =================================================================
    # Ion Drift Time Table
    # =================================================================
    println("\n--- Ion Drift Times (v_ion = $(v_drift_default) cm/s) ---")
    println("\nL (cm)  |  t_drift (s)")
    println("-" ^ 25)

    for L_cm in [50.0, 75.0, 100.0, 125.0, 150.0]
        t_s = ion_drift_time_s(L_cm)
        @printf("  %3.0f   |    %.2f\n", L_cm, t_s)
    end

    println("\nDone!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
