#!/usr/bin/env julia

"""
MARS Final Design — Complete Module
=====================================
Longerons: 3 × 5 mm OD, 0.5 mm wall, Ti Gr.5
Skin: Kapton, 50 µm
Belt: Kapton, 125 µm × 10 mm wide
Cable: Kapton flat cable
Ribs: 10 × HDPE formers
Ion plate: HDPE box 160×160×20 mm, 2 mm wall
Motor torque: 40 N·m
"""

# Activate project environment
const pdir = joinpath(dirname(@__FILE__), "..")
using Pkg
Pkg.activate(pdir)
Pkg.instantiate()

# ─── Constants ────────────────────────────────────────────────────────────

const ρ_Ti     = 4430.0    # kg/m³
const ρ_HDPE   = 960.0     # kg/m³
const ρ_Kapton = 1420.0    # kg/m³
const E_Ti     = 114e9     # Pa
const g        = 9.81      # m/s²

# ─── Tube geometry ────────────────────────────────────────────────────────

A_tube(d_o, t_w) = (π/4) * (d_o^2 - (d_o - 2t_w)^2)
I_tube(d_o, t_w) = (π/64) * (d_o^4 - (d_o - 2t_w)^4)

# ─── Linear mass density — component by component ────────────────────────

"""
    μ_longerons(N, d_o, t_w)

Linear mass density of N Ti Gr.5 longerons (kg/m).
"""
μ_longerons(N, d_o, t_w) = ρ_Ti * N * A_tube(d_o, t_w)

"""
    μ_skin(perimeter, t_skin)

Linear mass density of Kapton skin (kg/m).
    perimeter : NACA 0012 perimeter (m)
    t_skin    : skin thickness (m)
"""
μ_skin(perimeter, t_skin) = ρ_Kapton * perimeter * t_skin

"""
    μ_belt(width, thickness)

Linear mass density of Kapton belt (kg/m).
Belt runs the full arm length (×2 for loop).
"""
μ_belt(width, thickness) = ρ_Kapton * width * thickness * 2

"""
    μ_cable(width, thickness)

Linear mass density of Kapton flat cable (kg/m).
"""
μ_cable(width, thickness) = ρ_Kapton * width * thickness

"""
    μ_ribs(n_ribs, m_rib, R)

Effective linear mass density of discrete HDPE rib formers (kg/m).
    n_ribs : number of ribs per arm
    m_rib  : mass per rib (kg)
    R      : arm length (m)
"""
μ_ribs(n_ribs, m_rib, R) = n_ribs * m_rib / R

"""
    μ_total(; N, d_o, t_w, perimeter, t_skin, belt_w, belt_t,
              cable_w, cable_t, n_ribs, m_rib, R)

Total effective linear mass density of one arm (kg/m).
"""
function μ_total(; N=3, d_o=5e-3, t_w=0.5e-3,
                   perimeter=0.165, t_skin=50e-6,
                   belt_w=10e-3, belt_t=125e-6,
                   cable_w=10e-3, cable_t=50e-6,
                   n_ribs=10, m_rib=5e-3, R=1.6)
    μ_l = μ_longerons(N, d_o, t_w)
    μ_s = μ_skin(perimeter, t_skin)
    μ_b = μ_belt(belt_w, belt_t)
    μ_c = μ_cable(cable_w, cable_t)
    μ_r = μ_ribs(n_ribs, m_rib, R)
    μ_l + μ_s + μ_b + μ_c + μ_r
end

# ─── Ion plate mass (HDPE box) ────────────────────────────────────────────

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

function wake_velocity(R; d_wake=9.6e-3, Δθ=π/2, τ_motor=40.0, μ=0.14, m_plate=0.15)
    t_rot = rotation_time(μ, m_plate, R; τ_motor)
    α = 4Δθ / t_rot^2
    sqrt(α * R * d_wake)
end

function settling_time(R; v_d=0.10, n=1.2, d_wake=9.6e-3,
                       Δθ=π/2, τ_motor=40.0, μ=0.14, m_plate=0.15)
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

# ─── Bending and deflection ───────────────────────────────────────────────

M_tang(μ, α, R) = μ * α * R^3 / (9√3)
M_grav(μ, R) = μ * g * R^2 / 8
σ_bending(M, Z) = M / Z
δ_tang_max(μ, α, R, EI) = μ * α * R^5 / (153 * EI)
δ_grav_max(μ, R, EI) = 5 * μ * g * R^4 / (384 * EI)

# ─── Full summary ─────────────────────────────────────────────────────────

function mars_full_summary()
    # Design parameters
    N = 3; d_o = 5e-3; t_w = 0.5e-3; R = 1.6; L_fid = 1.50
    pos = [-30e-3, 0.0, 30e-3]
    τ_motor = 40.0

    # ── Arm mass budget ──
    μ_l = μ_longerons(N, d_o, t_w)
    μ_s = μ_skin(0.165, 50e-6)
    μ_b = μ_belt(10e-3, 125e-6)
    μ_c = μ_cable(10e-3, 50e-6)
    μ_r = μ_ribs(10, 5e-3, R)
    μ  = μ_l + μ_s + μ_b + μ_c + μ_r
    ma = μ * R

    println("=" ^ 60)
    println("  MARS Final Design Summary")
    println("=" ^ 60)
    println()
    println("── Arm mass budget (per arm) ──")
    println("  Ti longerons (3×5/0.5)  : μ = ", round(μ_l, digits=4), " kg/m  → ", round(μ_l*R, digits=3), " kg")
    println("  Kapton skin (50 µm)     : μ = ", round(μ_s, digits=4), " kg/m  → ", round(μ_s*R, digits=3), " kg")
    println("  Kapton belt (125 µm)    : μ = ", round(μ_b, digits=4), " kg/m  → ", round(μ_b*R, digits=3), " kg")
    println("  Kapton flat cable (50µm): μ = ", round(μ_c, digits=4), " kg/m  → ", round(μ_c*R, digits=3), " kg")
    println("  HDPE ribs (10×5 g)      : μ = ", round(μ_r, digits=4), " kg/m  → ", round(μ_r*R, digits=3), " kg")
    println("  ────────────────────────────────────────────")
    println("  Total                   : μ = ", round(μ, digits=4), " kg/m  → ", round(ma, digits=3), " kg")
    println()

    # ── Ion plate ──
    mp = m_plate_box(0.16, 0.16, 0.02, 0.002)
    println("── Ion plate (HDPE box 160×160×20, 2 mm wall + 30 g CMOS/screws) ──")
    println("  m_plate    = ", round(mp, digits=3), " kg")
    println()

    # ── Section properties (longerons only for structural calc) ──
    Is = I_section(N, d_o, t_w, pos)
    Zs = Z_section(N, d_o, t_w, pos)
    EI = EI_flex(N, d_o, t_w, pos)

    # ── Inertia (use total μ for dynamics) ──
    Ia = I_arms(μ, R)
    Ip = I_plates(mp, R)
    It = Ia + Ip
    II = inertia_term(μ, mp, R)
    DD = drag_term(R)

    println("── Moments of inertia ──")
    println("  I_arms     = ", round(Ia, digits=3), " kg·m²")
    println("  I_plates   = ", round(Ip, digits=3), " kg·m²")
    println("  I_total    = ", round(It, digits=3), " kg·m²")
    println("  𝓘          = ", round(II, digits=2), " N·m·s²")
    println("  𝓓          = ", round(DD, digits=2), " N·m·s²")
    println("  𝓓/𝓘        = ", round(100DD/II, digits=1), "%")
    println()

    # ── Rotation ──
    t_rot = rotation_time(μ, mp, R; τ_motor)
    α = 4(π/2) / t_rot^2
    ω_max = 2(π/2) / t_rot
    v_tip = ω_max * R

    println("── Rotation (τ_motor = ", τ_motor, " N·m) ──")
    println("  t_rot      = ", round(t_rot, digits=2), " s")
    println("  α          = ", round(α, digits=1), " rad/s²")
    println("  ω_max      = ", round(ω_max, digits=1), " rad/s")
    println("  v_tip      = ", round(v_tip, digits=1), " m/s")
    println()

    # ── Settling ──
    u₀ = wake_velocity(R; τ_motor, μ, m_plate=mp)
    τ₀ = 9.6e-3 / u₀
    t_set = settling_time(R; τ_motor=τ_motor, μ=μ, m_plate=mp)

    println("── Settling ──")
    println("  u₀         = ", round(u₀, digits=2), " m/s")
    println("  τ₀         = ", round(τ₀*1e3, digits=1), " ms")
    println("  t_settle   = ", round(t_set, digits=2), " s")
    println()

    # ── Performance ──
    t_cycle = t_rot + t_set
    Zd = 0.10 * t_cycle
    ε = (L_fid - Zd) / L_fid

    println("── Performance ──")
    println("  t_cycle    = ", round(t_cycle, digits=2), " s")
    println("  Z_dead     = ", round(Zd*100, digits=1), " cm")
    println("  ε_geo      = ", round(100ε, digits=1), "%")
    println()

    # ── Structural (longerons only carry load) ──
    Mt = M_tang(μ, α, R)
    Mg = M_grav(μ, R)
    σt = σ_bending(Mt, Zs)
    σg = σ_bending(Mg, Zs)
    σ_allow = 400e6
    dt = abs(δ_tang_max(μ, α, R, EI))
    dg = δ_grav_max(μ, R, EI)

    println("── Structural ──")
    println("  I_section  = ", round(Is*1e12, digits=0), " mm⁴")
    println("  Z_section  = ", round(Zs*1e9, digits=0), " mm³")
    println("  EI         = ", round(EI, digits=0), " N·m²")
    println("  M_tang     = ", round(Mt, digits=2), " N·m")
    println("  M_grav     = ", round(Mg, digits=2), " N·m")
    println("  σ_total    = ", round((σt+σg)/1e6, digits=1), " MPa")
    println("  σ_allow    = ", round(σ_allow/1e6, digits=0), " MPa")
    println("  margin     = ", round(σ_allow/(σt+σg), digits=0), "×")
    println("  δ_tang     = ", round(dt*1e3, digits=2), " mm")
    println("  δ_grav     = ", round(dg*1e3, digits=2), " mm")
    println("=" ^ 60)
end

if abspath(PROGRAM_FILE) == @__FILE__
    mars_full_summary()
end