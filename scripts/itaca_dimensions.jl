#!/usr/bin/env julia

"""
ITACA Full Detector Geometry
============================
Computes and prints complete detector dimensions including radial and axial structure.
"""

include("itaca_f.jl")

function main()
    # Compute full detector dimensions with default parameters
    dims = itaca_full_dimensions(; D_fid_m=D_fid_default, L_TPC_m=L_fid_default)

    # Print formatted summary
    print_itaca_dimensions(dims)

    # Additional summary
    println("\n--- Summary ---")
    @printf("  Pressure vessel: %.2f m (ID) × %.2f m (H)\n", dims.D_PV_inner_m, dims.L_internal_m)
    @printf("  Wall thickness:  %.1f mm (cylinder), %.1f mm (hemisphere)\n", dims.t_cyl_mm, dims.t_hem_mm)

    # Fiducial volume and mass
    V_fid = π * (dims.D_fid_m/2)^2 * dims.L_TPC_m
    m_fid = V_fid * ρ_Xe_15bar
    @printf("  Fiducial volume: %.1f m³\n", V_fid)
    @printf("  Fiducial mass:   %.0f kg (%.2f ton)\n", m_fid, m_fid/1000)

    println("\nDone!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
