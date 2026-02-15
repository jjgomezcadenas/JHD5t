#!/usr/bin/env julia

"""
MARS ITACA — Unified Kinematics & Structural Analysis
======================================================
Longerons: 3 × 5 mm OD, 0.5 mm wall, Ti Gr.5
Ion plate: HDPE box 160×160×20 mm, 2 mm wall
Motor torque: 40 N·m (reduced to limit wake)
"""

# Activate project environment
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

using Printf

# ─── Material properties ──────────────────────────────────────────────────

const ρ_Ti   = 4430.0    # kg/m³
const ρ_HDPE = 960.0     # kg/m³
const E_Ti   = 114e9     # Pa
const g      = 9.81      # m/s²

# ─── Tube geometry ────────────────────────────────────────────────────────

A_tube(d_o, t_w) = (π/4) * (d_o^2 - (d_o - 2t_w)^2)
I_tube(d_o, t_w) = (π/64) * (d_o^4 - (d_o - 2t_w)^4)

# ─── Arm properties ───────────────────────────────────────────────────────

μ_arm(N, d_o, t_w; ρ=ρ_Ti) = ρ * N * A_tube(d_o, t_w)
m_arm(N, d_o, t_w, R; ρ=ρ_Ti) = μ_arm(N, d_o, t_w; ρ) * R

# ─── Ion plate mass (HDPE box) ────────────────────────────────────────────

"""
    m_plate_box(Lx, Ly, Lz, t_p; ρ=ρ_HDPE, m_extra=0.03)

Mass of an HDPE box structure (kg).
    Lx, Ly : footprint dimensions (m)
    Lz     : box height (m)
    t_p    : wall thickness (m)
    m_extra: CMOS + screws (kg)
"""
function m_plate_box(Lx, Ly, Lz, t_p; ρ=ρ_HDPE, m_extra=0.03)
    A_top_bot = 2 * Lx * Ly
    A_sides   = 2 * (Lx + Ly) * Lz
    V = (A_top_bot + A_sides) * t_p
    ρ * V + m_extra
end

# ─── Section properties ───────────────────────────────────────────────────

function I_section(N, d_o, t_w, positions)
    It = I_tube(d_o, t_w)
    At = A_tube(d_o, t_w)
    N * It + At * sum(y^2 for y in positions)
end

Z_section(N, d_o, t_w, positions) = I_section(N, d_o, t_w, positions) / maximum(abs.(positions))
EI_flex(N, d_o, t_w, positions; E=E_Ti) = E * I_section(N, d_o, t_w, positions)

# ─── Moments of inertia ──────────────────────────────────────────────────

I_arms(μ, R) = (2/3) * μ * R^3
I_plates(m_plate, R) = 2 * m_plate * R^2
I_total(μ, m_plate, R) = I_arms(μ, R) + I_plates(m_plate, R)

# ─── Torque balance ───────────────────────────────────────────────────────

inertia_term(μ, m_plate, R; Δθ=π/2) = 4 * I_total(μ, m_plate, R) * Δθ
drag_term(R; ρ=87.0, Cd=0.1, d_wake=9.6e-3, Δθ=π/2) = ρ * Cd * d_wake * R^4 * Δθ^2

function rotation_time(μ, m_plate, R;
                       τ_motor=40.0, ρ=87.0, Cd=0.1, d_wake=9.6e-3, Δθ=π/2)
    𝓘 = inertia_term(μ, m_plate, R; Δθ)
    𝓓 = drag_term(R; ρ, Cd, d_wake, Δθ)
    sqrt((𝓘 + 𝓓) / τ_motor)
end

# ─── Wake and settling ────────────────────────────────────────────────────

function wake_velocity(R; d_wake=9.6e-3, Δθ=π/2, τ_motor=40.0, μ=0.094, m_plate=0.15)
    t_rot = rotation_time(μ, m_plate, R; τ_motor)
    α = 4Δθ / t_rot^2
    sqrt(α * R * d_wake)
end

function settling_time(R; v_d=0.10, n=1.2, d_wake=9.6e-3,
                       Δθ=π/2, τ_motor=40.0, μ=0.094, m_plate=0.15)
    u₀ = wake_velocity(R; d_wake, Δθ, τ_motor, μ, m_plate)
    τ₀ = d_wake / u₀
    τ₀ * ((u₀ / v_d)^(2/n) - 1)
end

# ─── Dead zone and efficiency ─────────────────────────────────────────────

function dead_zone(μ, m_plate, R; v_d=0.10, n=1.2, d_wake=9.6e-3,
                   Δθ=π/2, τ_motor=40.0, ρ=87.0, Cd=0.1)
    t_rot = rotation_time(μ, m_plate, R; τ_motor, ρ, Cd, d_wake, Δθ)
    t_set = settling_time(R; v_d, n, d_wake, Δθ, τ_motor, μ, m_plate)
    v_d * (t_rot + t_set)
end

function fiducial_efficiency(μ, m_plate, R, L_fid; v_d=0.10, n=1.2, d_wake=9.6e-3,
                             Δθ=π/2, τ_motor=40.0, ρ=87.0, Cd=0.1)
    Z = dead_zone(μ, m_plate, R; v_d, n, d_wake, Δθ, τ_motor, ρ, Cd)
    (L_fid - Z) / L_fid
end

# ─── Bending moments ─────────────────────────────────────────────────────

M_tang(μ, α, R) = μ * α * R^3 / (9√3)
M_grav(μ, R) = μ * g * R^2 / 8
σ_bending(M, Z) = M / Z

# ─── Deflections ──────────────────────────────────────────────────────────

δ_tang_max(μ, α, R, EI) = μ * α * R^5 / (153 * EI)
δ_grav_max(μ, R, EI) = 5 * μ * g * R^4 / (384 * EI)

# ─── Full summary ─────────────────────────────────────────────────────────

function main()
    # Design parameters
    N = 3; d_o = 5e-3; t_w = 0.5e-3; R = 1.6; L_fid = 1.50
    pos = [-30e-3, 0.0, 30e-3]
    τ_motor = 40.0
    v_d = 0.10  # m/s (ion drift velocity)

    # Arm
    μ  = μ_arm(N, d_o, t_w)
    ma = m_arm(N, d_o, t_w, R)

    # Plate
    mp = m_plate_box(0.16, 0.16, 0.02, 0.002)

    # Section
    Is = I_section(N, d_o, t_w, pos)
    Zs = Z_section(N, d_o, t_w, pos)
    EI = EI_flex(N, d_o, t_w, pos)

    # Inertia
    Ia = I_arms(μ, R)
    Ip = I_plates(mp, R)
    It = Ia + Ip
    II = inertia_term(μ, mp, R)
    DD = drag_term(R)

    # Rotation
    t_rot = rotation_time(μ, mp, R; τ_motor)
    α = 4(π/2) / t_rot^2
    ω_max = 2(π/2) / t_rot
    v_tip = ω_max * R

    # Settling
    u₀ = wake_velocity(R; τ_motor, μ, m_plate=mp)
    τ₀ = 9.6e-3 / u₀
    t_set = settling_time(R; v_d, τ_motor, μ, m_plate=mp)

    # Performance
    t_cycle = t_rot + t_set
    Zd = dead_zone(μ, mp, R; v_d, τ_motor)
    ε = fiducial_efficiency(μ, mp, R, L_fid; v_d, τ_motor)

    # Structural
    Mt = M_tang(μ, α, R)
    Mg = M_grav(μ, R)
    σt = σ_bending(Mt, Zs)
    σg = σ_bending(Mg, Zs)
    σ_allow = 400e6
    dt = abs(δ_tang_max(μ, α, R, EI))
    dg = δ_grav_max(μ, R, EI)

    println("\n" * "=" ^ 60)
    println("          MARS ITACA — Full Design Summary")
    println("=" ^ 60)

    println("\n── Input Parameters ──")
    @printf("  Arm radius:           R = %.1f m\n", R)
    @printf("  Fiducial length:      L_fid = %.1f m (%.0f cm)\n", L_fid, L_fid * 100)
    @printf("  Ion drift velocity:   v_d = %.2f m/s (%.0f cm/s)\n", v_d, v_d * 100)
    @printf("  Motor torque:         τ_motor = %.0f N·m\n", τ_motor)

    println("\n── Arm Structure ──")
    @printf("  Longerons:            %d × %.1f mm OD, %.1f mm wall (Ti Gr.5)\n", N, d_o*1e3, t_w*1e3)
    @printf("  A_tube:               %.1f mm²\n", A_tube(d_o, t_w)*1e6)
    @printf("  μ (linear density):   %.3f kg/m\n", μ)
    @printf("  m_arm (per arm):      %.2f kg\n", ma)

    println("\n── Ion Plate ──")
    @printf("  HDPE box:             160×160×20 mm, 2 mm wall\n")
    @printf("  m_plate:              %.2f kg\n", mp)

    println("\n── Moment of Inertia ──")
    @printf("  I_arms = (2/3) μ R³:  %.4f kg·m²\n", Ia)
    @printf("  I_plates = 2 m R²:    %.3f kg·m²\n", Ip)
    @printf("  I_total:              %.3f kg·m²\n", It)
    @printf("  Inertia term 𝓘:       %.2f N·m·s²\n", II)
    @printf("  Drag term 𝓓:          %.2f N·m·s²  (%.1f%% of 𝓘)\n", DD, 100DD/II)

    println("\n── Rotation Dynamics ──")
    @printf("  t_rot:                %.2f s\n", t_rot)
    @printf("  α (acceleration):     %.1f rad/s²\n", α)
    @printf("  ω_max:                %.1f rad/s\n", ω_max)
    @printf("  v_tip:                %.1f m/s\n", v_tip)

    println("\n── Wake & Settling ──")
    @printf("  Wake velocity u₀:     %.2f m/s\n", u₀)
    @printf("  Ratio u₀/v_d:         %.1f\n", u₀ / v_d)
    @printf("  τ₀:                   %.1f ms\n", τ₀*1e3)
    @printf("  t_settle:             %.2f s\n", t_set)

    println("\n── MARS Performance ──")
    @printf("  Cycle time t_cycle:   %.2f s  (t_rot + t_settle)\n", t_cycle)
    @printf("  Dead zone Z_dead:     %.1f cm\n", Zd * 100)
    @printf("  Fiducial efficiency:  %.1f %%\n", ε * 100)

    println("\n── Structural Analysis ──")
    @printf("  I_section:            %.0f mm⁴\n", Is*1e12)
    @printf("  Z_section:            %.0f mm³\n", Zs*1e9)
    @printf("  EI (flexural):        %.0f N·m²\n", EI)
    @printf("  M_tang (tangential):  %.2f N·m\n", Mt)
    @printf("  M_grav (gravity):     %.2f N·m\n", Mg)
    @printf("  σ_total:              %.1f MPa\n", (σt+σg)/1e6)
    @printf("  σ_allow:              %.0f MPa\n", σ_allow/1e6)
    @printf("  Safety margin:        %.0f×\n", σ_allow/(σt+σg))
    @printf("  δ_tang:               %.2f mm\n", dt*1e3)
    @printf("  δ_grav:               %.2f mm\n", dg*1e3)

    # Sensitivity table
    println("\n── Sensitivity to v_d ──")
    println("\n  v_d (m/s)  |  t_settle (s)  |  Z_dead (cm)  |  ε_geo (%)")
    println("  " * "-" ^ 55)
    for v in [0.10, 0.12, 0.15, 0.18, 0.20]
        t_s = settling_time(R; v_d=v, τ_motor, μ, m_plate=mp)
        Z = dead_zone(μ, mp, R; v_d=v, τ_motor)
        ε_v = fiducial_efficiency(μ, mp, R, L_fid; v_d=v, τ_motor)
        @printf("     %.2f     |      %.2f       |     %.1f      |    %.1f\n", v, t_s, Z * 100, ε_v * 100)
    end

    println("\n" * "=" ^ 60)
    println("Done!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
