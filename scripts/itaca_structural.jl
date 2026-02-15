#!/usr/bin/env julia

"""
ITACA MARS Arm Structural Analysis
==================================
Structural verification of MARS arm: geometry, bending, stress, deflection.
"""

include("itaca_f.jl")

function main()
    println("\n" * "=" ^ 60)
    println("          MARS Arm Structural Analysis")
    println("=" ^ 60)

    # Arm parameters: 3 tubes, 8 mm OD, 0.8 mm wall, at y = -30, 0, +30 mm
    N = 3
    d_o = 8e-3       # m (outer diameter)
    t_w = 0.8e-3     # m (wall thickness)
    R = 1.6          # m (arm radius)
    pos = [-30e-3, 0.0, 30e-3]  # m (longeron positions)

    m_plate = 0.6    # kg (ion plate mass)

    # =================================================================
    # Tube and Arm Geometry
    # =================================================================
    println("\n--- Tube Geometry ---")
    @printf("  Number of longerons:  N = %d\n", N)
    @printf("  Outer diameter:       d_o = %.1f mm\n", d_o * 1e3)
    @printf("  Wall thickness:       t_w = %.1f mm\n", t_w * 1e3)
    @printf("  Inner diameter:       d_i = %.1f mm\n", (d_o - 2t_w) * 1e3)

    At = A_tube(d_o, t_w)
    It = I_tube(d_o, t_w)
    @printf("  Tube area:            A_tube = %.1f mm²\n", At * 1e6)
    @printf("  Tube moment:          I_tube = %.1f mm⁴\n", It * 1e12)

    println("\n--- Arm Parameters ---")
    @printf("  Arm radius:           R = %.1f m\n", R)
    @printf("  Longeron positions:   y = [%.0f, %.0f, %.0f] mm\n", pos[1]*1e3, pos[2]*1e3, pos[3]*1e3)

    μ = μ_arm(N, d_o, t_w)
    m = m_arm(N, d_o, t_w, R)
    Is = I_section(N, d_o, t_w, pos)
    Zs = Z_section(N, d_o, t_w, pos)

    @printf("  Linear mass density:  μ = %.3f kg/m\n", μ)
    @printf("  Arm mass:             m_arm = %.2f kg\n", m)
    @printf("  Section moment:       I_section = %.0f mm⁴\n", Is * 1e12)
    @printf("  Section modulus:      Z_section = %.0f mm³\n", Zs * 1e9)

    # =================================================================
    # Load Case Analysis
    # =================================================================
    # Angular acceleration from MARS kinematics
    t_rot = rotation_time(μ, m_plate, R)
    α = 4 * (π/2) / t_rot^2

    println("\n--- Loading Conditions ---")
    @printf("  Rotation time:        t_rot = %.2f s\n", t_rot)
    @printf("  Angular acceleration: α = %.1f rad/s²\n", α)
    @printf("  Gravity:              g = %.2f m/s²\n", g)

    # =================================================================
    # Bending Moments
    # =================================================================
    Mt = M_tang(μ, α, R)
    Mg = M_grav(μ, R)
    Mtot = Mt + Mg

    println("\n--- Bending Moments ---")
    @printf("  Tangential (inertia): M_tang = %.2f N·m  (at r = R/√3 = %.2f m)\n", Mt, r_tang_max(R))
    @printf("  Gravity (uniform):    M_grav = %.2f N·m  (at r = R/2 = %.2f m)\n", Mg, R/2)
    @printf("  Total (worst case):   M_total = %.2f N·m\n", Mtot)

    # =================================================================
    # Stress Analysis
    # =================================================================
    σt = σ_bending(Mt, Zs)
    σg = σ_bending(Mg, Zs)
    σtot = σ_bending(Mtot, Zs)

    σ_yield = 830e6   # Pa, Ti Gr.5 yield strength
    σ_allow = 400e6   # Pa, with safety factor ~2

    println("\n--- Stress Analysis ---")
    @printf("  σ_tang:               %.1f MPa\n", σt / 1e6)
    @printf("  σ_grav:               %.1f MPa\n", σg / 1e6)
    @printf("  σ_total:              %.1f MPa\n", σtot / 1e6)
    @printf("  σ_yield (Ti Gr.5):    %.0f MPa\n", σ_yield / 1e6)
    @printf("  σ_allow (SF=2):       %.0f MPa\n", σ_allow / 1e6)
    @printf("  Safety margin:        %.0f×\n", σ_allow / σtot)

    # =================================================================
    # Deflection Analysis
    # =================================================================
    EI = EI_flex(N, d_o, t_w, pos)
    dt = abs(δ_tang_max(μ, α, R, EI))
    dg = δ_grav_max(μ, R, EI)
    dtot = dt + dg

    println("\n--- Deflection Analysis ---")
    @printf("  Flexural rigidity:    EI = %.0f N·m²\n", EI)
    @printf("  Young's modulus:      E = %.0f GPa (Ti Gr.5)\n", E_Ti / 1e9)
    @printf("  δ_tang (transient):   %.2f mm  (at r ≈ 0.52R)\n", dt * 1e3)
    @printf("  δ_grav (static):      %.2f mm  (at r = R/2)\n", dg * 1e3)
    @printf("  δ_total (worst):      %.2f mm\n", dtot * 1e3)

    # Deflection limit (typically L/300 for structural members)
    δ_limit = R / 300
    @printf("  Deflection limit:     %.1f mm  (R/300)\n", δ_limit * 1e3)
    @printf("  Margin:               %.0f×\n", δ_limit / dtot)

    println("\nDone!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
