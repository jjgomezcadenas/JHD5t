#!/usr/bin/env julia

"""
ITACA Dense Silicon Plane (DSP) Analysis
=========================================
Computes DSP parameters: number of SiPM channels, coverage, and areas.
"""

include("itaca_f.jl")

function main()
    println("\n" * "=" ^ 60)
    println("          ITACA Dense Silicon Plane (DSP)")
    println("=" ^ 60)

    # DSP diameter = TPC fiducial diameter
    D_cm = D_fid_default * 100.0  # m to cm

    # Print full summary
    print_dsp_summary(D_cm)

    # Additional derived values
    n_channels = dsp_n_channels(D_cm)
    coverage = dsp_coverage()
    total_sipm_area = dsp_total_sipm_area_m2(D_cm)
    dsp_area = π * (D_cm / 200.0)^2  # m²

    println("\n--- Derived Parameters ---")
    @printf("  Channels per row:   ~%d (along diameter)\n", dsp_n_sipms_1d(D_cm))
    @printf("  SiPM area per ch:   %.0f mm²\n", SiPM_size_mm^2)
    @printf("  Pitch area:         %.0f mm²\n", DSP_pitch_mm^2)

    # Estimate power (typical values)
    power_per_sipm_mW = 5.0  # mW per channel (typical)
    total_power_W = n_channels * power_per_sipm_mW / 1e3

    println("\n--- System Estimates ---")
    @printf("  Power (@ 5 mW/ch):  %.0f W\n", total_power_W)

    println("\n" * "-" ^ 60)
    println("Done!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
