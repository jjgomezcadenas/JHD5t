#!/usr/bin/env julia

"""
itaca Detector Script

This script creates the itaca detector configuration for background calculations.
The detector includes copper shielding, PTFE lining, and a central cathode grid.

Usage:
    julia itaca.jl [options]

Options:
    --plots         Create and save detector plots (default: false)
    --dataframe     Save activity summary to CSV (default: true)
    --no-dataframe  Disable saving activity summary to CSV
    --print=MODE    Print mode: "concise" or "verbose" (default: concise)
    --help, -h      Show this help message
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JHD5t
using Unitful
using Printf
using Plots
using DataFrames
using CSV
using ArgParse

# Include MARS parameters (materials, masses, kinematics)
include("mars_pars.jl")

# Fiducial dimensions
const L_fid_default = 1.5   # m (fiducial length)
const D_fid_default = 3.2   # m (fiducial diameter)
const BFD_thickness_cm = 1.0        # Buffer Field region
const clearance_cm = 1.0            # Clearances
const ICS_thickness_cm = 3.0       # Inner Copper Shield

# Axial components
const DSP_LG_FATGEM_cm = 10.0       # DSP + Light Guide + FAT-GEM combined
const MRS_thickness_cm = 10.0       # MRS area (below cathode)

# ─── MARS Component Mass Functions ───────────────────────────────────────────
# These functions return total mass (kg) for MARS components using parameters
# from mars_pars.jl. Default values match the MARS final design.

"""
    mass_longerons(; N=3, d_o=5e-3, t_w=0.5e-3, R=1.6)

Total mass of Ti Gr.5 longerons for both arms (kg).
Default: 3 longerons, 5 mm OD, 0.5 mm wall, 1.6 m arm length.
"""
function mass_longerons(; N=3, d_o=5e-3, t_w=0.5e-3, R=1.6)
    μ = μ_longerons(N, d_o, t_w)
    return 2 * μ * R  # Two arms
end

"""
    mass_mars_skin(; perimeter=0.165, t_skin=50e-6, R=1.6)

Total mass of Kapton skin for both arms (kg).
Default: NACA 0012 perimeter 165 mm, 50 µm thickness, 1.6 m arm length.
"""
function mass_mars_skin(; perimeter=0.165, t_skin=50e-6, R=1.6)
    μ = μ_skin(perimeter, t_skin)
    return 2 * μ * R  # Two arms
end

"""
    mass_mars_belt(; width=10e-3, thickness=125e-6, R=1.6)

Total mass of Kapton belt for both arms (kg).
Default: 10 mm wide, 125 µm thick, 1.6 m arm length.
Note: Belt runs full length ×2 for loop (accounted in μ_belt).
"""
function mass_mars_belt(; width=10e-3, thickness=125e-6, R=1.6)
    μ = μ_belt(width, thickness)
    return 2 * μ * R  # Two arms
end

"""
    mass_mars_cables(; width=10e-3, thickness=50e-6, R=1.6)

Total mass of Kapton flat cables for both arms (kg).
Default: 10 mm wide, 50 µm thick, 1.6 m arm length.
"""
function mass_mars_cables(; width=10e-3, thickness=50e-6, R=1.6)
    μ = μ_cable(width, thickness)
    return 2 * μ * R  # Two arms
end

"""
    mass_ion_carrier(; Lx=0.16, Ly=0.16, Lz=0.02, t_p=0.002, m_extra=0.03)

Total mass of both ion carrier plates (kg).
Default: HDPE box 160×160×20 mm, 2 mm wall, 30 g CMOS/screws.
"""
function mass_ion_carrier(; Lx=0.16, Ly=0.16, Lz=0.02, t_p=0.002, m_extra=0.03)
    return 2 * m_plate_box(Lx, Ly, Lz, t_p; m_extra=m_extra)  # Two plates
end

"""
    mass_mars_ribs(; n_ribs=10, m_rib=5e-3)

Total mass of HDPE rib formers for both arms (kg).
Default: 10 ribs per arm, 5 g each.
"""
function mass_mars_ribs(; n_ribs=10, m_rib=5e-3)
    return 2 * n_ribs * m_rib  # Two arms
end

"""
    mass_mars_total(; N=3, d_o=5e-3, t_w=0.5e-3, R=1.6,
                      perimeter=0.165, t_skin=50e-6,
                      belt_w=10e-3, belt_t=125e-6,
                      cable_w=10e-3, cable_t=50e-6,
                      n_ribs=10, m_rib=5e-3,
                      Lx=0.16, Ly=0.16, Lz=0.02, t_p=0.002, m_extra=0.03)

Total mass of the complete MARS system (kg).
Includes: longerons, skin, belt, cables, ribs, and ion carriers for both arms.
"""
function mass_mars_total(; N=3, d_o=5e-3, t_w=0.5e-3, R=1.6,
                           perimeter=0.165, t_skin=50e-6,
                           belt_w=10e-3, belt_t=125e-6,
                           cable_w=10e-3, cable_t=50e-6,
                           n_ribs=10, m_rib=5e-3,
                           Lx=0.16, Ly=0.16, Lz=0.02, t_p=0.002, m_extra=0.03)
    m_long = mass_longerons(; N, d_o, t_w, R)
    m_skin = mass_mars_skin(; perimeter, t_skin, R)
    m_belt = mass_mars_belt(; width=belt_w, thickness=belt_t, R)
    m_cables = mass_mars_cables(; width=cable_w, thickness=cable_t, R)
    m_ribs = mass_mars_ribs(; n_ribs, m_rib)
    m_plates = mass_ion_carrier(; Lx, Ly, Lz, t_p, m_extra)
    return m_long + m_skin + m_belt + m_cables + m_ribs + m_plates
end

# ─── DSP Honeycomb Mass Functions ───────────────────────────────────────────

"""
    honeycomb_fill_factor(d_channel, pitch)

Calculate the fill factor for a hexagonal honeycomb structure with circular channels.

For a hexagonal lattice with pitch p (center-to-center distance) and channel
diameter d, the unit cell is a hexagon with area A_hex = (√3/2) × p².
Each channel removes area A_hole = π × (d/2)².

Fill factor f = 1 - (A_hole / A_hex) = 1 - (π/2√3) × (d/p)²

# Arguments
- `d_channel`: Channel diameter (m)
- `pitch`: Center-to-center distance between channels (m)

# Returns
- Fill factor (dimensionless, between 0 and 1)
"""
function honeycomb_fill_factor(d_channel, pitch)
    ratio = d_channel / pitch
    return 1.0 - (π / (2 * √3)) * ratio^2
end

"""
    mass_dsp_honeycomb(; D=2.30, L=0.10, d_channel=6e-3, pitch=10e-3, ρ=ρ_HDPE)

Calculate the effective mass of the DSP honeycomb structure.

The DSP is a disk with a honeycomb structure made of circular channels at
a hexagonal pitch. The mass is computed as:
    mass = (solid disk volume) × fill_factor × density

# Arguments
- `D`: Disk diameter (m), default 2.30 m (230 cm)
- `L`: Disk thickness (m), default 0.10 m (10 cm)
- `d_channel`: Channel diameter (m), default 6e-3 (6 mm)
- `pitch`: Channel pitch center-to-center (m), default 10e-3 (10 mm)
- `ρ`: Material density (kg/m³), default ρ_HDPE (960 kg/m³)

# Returns
- Named tuple: (mass_solid, mass_honeycomb, fill_factor, volume_solid, volume_effective)
  - mass_solid: Mass if disk were solid (kg)
  - mass_honeycomb: Effective mass with honeycomb (kg)
  - fill_factor: Fill factor (dimensionless)
  - volume_solid: Solid disk volume (m³)
  - volume_effective: Effective volume with honeycomb (m³)
"""
function mass_dsp_honeycomb(; D=2.30, L=0.10, d_channel=6e-3, pitch=10e-3, ρ=ρ_HDPE)
    # Solid disk volume
    R = D / 2.0
    V_solid = π * R^2 * L

    # Fill factor
    f = honeycomb_fill_factor(d_channel, pitch)

    # Effective volume and mass
    V_eff = V_solid * f
    m_solid = V_solid * ρ
    m_honeycomb = V_eff * ρ

    return (
        mass_solid = m_solid,
        mass_honeycomb = m_honeycomb,
        fill_factor = f,
        volume_solid = V_solid,
        volume_effective = V_eff
    )
end

"""
    print_dsp_honeycomb_summary(; D=2.30, L=0.10, d_channel=6e-3, pitch=10e-3, ρ=ρ_HDPE)

Print a summary of the DSP honeycomb structure parameters and computed mass.
"""
function print_dsp_honeycomb_summary(; D=2.30, L=0.10, d_channel=6e-3, pitch=10e-3, ρ=ρ_HDPE)
    result = mass_dsp_honeycomb(; D, L, d_channel, pitch, ρ)

    println("\n" * "="^60)
    println("          DSP Honeycomb Structure Summary")
    println("="^60)

    println("\n── Geometry ──")
    @printf("  Disk diameter:     D = %.0f cm\n", D * 100)
    @printf("  Disk thickness:    L = %.0f cm\n", L * 100)
    @printf("  Channel diameter:  d = %.1f mm\n", d_channel * 1000)
    @printf("  Channel pitch:     p = %.1f mm\n", pitch * 1000)
    @printf("  d/p ratio:         %.2f\n", d_channel / pitch)

    println("\n── Material ──")
    @printf("  Density:           ρ = %.0f kg/m³\n", ρ)

    println("\n── Computed Values ──")
    @printf("  Fill factor:       f = %.3f (%.1f%% material)\n", result.fill_factor, result.fill_factor * 100)
    @printf("  Solid volume:      V = %.4f m³ = %.2f L\n", result.volume_solid, result.volume_solid * 1000)
    @printf("  Effective volume:  V_eff = %.4f m³ = %.2f L\n", result.volume_effective, result.volume_effective * 1000)
    @printf("  Mass if solid:     m_solid = %.2f kg\n", result.mass_solid)
    @printf("  Honeycomb mass:    m_hc = %.2f kg\n", result.mass_honeycomb)
    @printf("  Mass reduction:    %.1f%%\n", (1 - result.fill_factor) * 100)

    println("="^60)

    return result
end

# ─── Cathode Mesh Mass Functions ───────────────────────────────────────────

# Cathode mesh default parameters (from itaca_f.jl)
const w_mesh_default = 200e-6  # m (wire diameter, 200 μm)
const M_mesh_default = 5e-3    # m (mesh pitch, 5 mm)
const ρ_Fe316Ti = 7870.0       # kg/m³ (stainless steel 316Ti density)

"""
    mesh_fill_factor(w, M)

Calculate the fill factor for a square wire mesh.

For a square mesh with wire diameter w and pitch M:
- Each unit cell (M × M) contains wire of total length ≈ 2M
- Wire cross-section: A = π(w/2)²
- Fill factor = (wire volume per cell) / (cell volume) ≈ (π/2)(w/M)²

# Arguments
- `w`: Wire diameter (m)
- `M`: Mesh pitch (m)

# Returns
- Fill factor (dimensionless)
"""
function mesh_fill_factor(w, M)
    return (π / 2) * (w / M)^2
end

"""
    mass_cathode_mesh(; D=3.20, w=w_mesh_default, M=M_mesh_default, ρ=ρ_Fe316Ti)

Calculate the mass of the cathode wire mesh.

The cathode is a square wire mesh (Fe316Ti) covering a circular disk.
For a square mesh, wire length per unit area ≈ 2/M (two perpendicular sets).

# Arguments
- `D`: Disk diameter (m), default 3.20 m (320 cm, full TPC diameter)
- `w`: Wire diameter (m), default 200 μm
- `M`: Mesh pitch (m), default 5 mm
- `ρ`: Wire material density (kg/m³), default ρ_Fe316Ti (7870 kg/m³)

# Returns
- Named tuple: (mass, fill_factor, wire_length, disk_area, wire_diameter, pitch)
"""
function mass_cathode_mesh(; D=3.20, w=w_mesh_default, M=M_mesh_default, ρ=ρ_Fe316Ti)
    # Disk area
    R = D / 2.0
    A_disk = π * R^2  # m²

    # Wire cross-section
    A_wire = π * (w / 2)^2  # m²

    # Wire length per unit area (two perpendicular sets of wires)
    wire_per_m2 = 2.0 / M  # m/m²

    # Total wire length
    L_wire = wire_per_m2 * A_disk  # m

    # Total mass
    m = ρ * A_wire * L_wire  # kg

    # Fill factor
    f = mesh_fill_factor(w, M)

    return (
        mass = m,
        fill_factor = f,
        wire_length = L_wire,
        disk_area = A_disk,
        wire_diameter = w,
        pitch = M
    )
end

"""
    print_cathode_mesh_summary(; D=3.20, w=w_mesh_default, M=M_mesh_default, ρ=ρ_Fe316Ti)

Print a summary of the cathode mesh parameters and computed mass.
"""
function print_cathode_mesh_summary(; D=3.20, w=w_mesh_default, M=M_mesh_default, ρ=ρ_Fe316Ti)
    result = mass_cathode_mesh(; D, w, M, ρ)

    println("\n" * "="^60)
    println("          Cathode Mesh Summary (Fe316Ti)")
    println("="^60)

    println("\n── Geometry ──")
    @printf("  Disk diameter:     D = %.0f cm\n", D * 100)
    @printf("  Wire diameter:     w = %.0f μm\n", w * 1e6)
    @printf("  Mesh pitch:        M = %.1f mm\n", M * 1000)
    @printf("  w/M ratio:         %.4f\n", w / M)

    println("\n── Material ──")
    @printf("  Material:          Fe316Ti (Stainless Steel)\n")
    @printf("  Density:           ρ = %.0f kg/m³\n", ρ)

    println("\n── Computed Values ──")
    @printf("  Fill factor:       f = %.6f (%.4f%%)\n", result.fill_factor, result.fill_factor * 100)
    @printf("  Disk area:         A = %.4f m²\n", result.disk_area)
    @printf("  Total wire length: L = %.1f m (%.1f km)\n", result.wire_length, result.wire_length / 1000)
    @printf("  Cathode mass:      m = %.3f kg (%.1f g)\n", result.mass, result.mass * 1000)

    println("="^60)

    return result
end

# ─── Activity Calculation Functions ────────────────────────────────────────

# Specific activities (Bq/kg) - from MaterialsData.jl via JHD5t materials module
# These are duplicated here for standalone use; the JHD5t module exports them with units
const SA_Bi214 = Dict(
    "Fe316Ti" => 1.9e-3,   # Bq/kg
    "Cu"      => 12.0e-6,  # Bq/kg
    "PTFE"    => 25.0e-6,  # Bq/kg
    "HDPE"    => 62.0e-6,  # Bq/kg (Poly)
    "Kapton"  => 1.0e-3,   # Bq/kg
    "Ti"      => 0.93e-3,  # Bq/kg
)

const SA_Tl208 = Dict(
    "Fe316Ti" => 0.4e-3,   # Bq/kg
    "Cu"      => 1.4e-6,   # Bq/kg
    "PTFE"    => 10.0e-6,  # Bq/kg
    "HDPE"    => 8.0e-6,   # Bq/kg (Poly)
    "Kapton"  => 1.0e-3,   # Bq/kg
    "Ti"      => 0.22e-3,  # Bq/kg
)

"""
    activity_bi214(mass_kg, material::String)

Compute Bi-214 activity (Bq) for a given mass and material.

# Arguments
- `mass_kg`: Mass in kg
- `material`: Material name (Fe316Ti, Cu, PTFE, HDPE, Kapton, Ti)

# Returns
- Activity in Bq
"""
function activity_bi214(mass_kg, material::String)
    if !haskey(SA_Bi214, material)
        error("Unknown material: $material. Valid: $(keys(SA_Bi214))")
    end
    return mass_kg * SA_Bi214[material]
end

"""
    activity_tl208(mass_kg, material::String)

Compute Tl-208 activity (Bq) for a given mass and material.

# Arguments
- `mass_kg`: Mass in kg
- `material`: Material name (Fe316Ti, Cu, PTFE, HDPE, Kapton, Ti)

# Returns
- Activity in Bq
"""
function activity_tl208(mass_kg, material::String)
    if !haskey(SA_Tl208, material)
        error("Unknown material: $material. Valid: $(keys(SA_Tl208))")
    end
    return mass_kg * SA_Tl208[material]
end

# ─── Per-Component Activity Functions ──────────────────────────────────────

"""
    activity_bfd(; R_in=1.60, R_out=1.61, L=1.70)

Compute activity of the BFD (Barrel Fiber Detector) - PTFE shell.

# Arguments
- `R_in`: Inner radius (m), default 1.60 m
- `R_out`: Outer radius (m), default 1.61 m
- `L`: Length (m), default 1.70 m

# Returns
- Named tuple: (mass_kg, bi214_Bq, tl208_Bq, bi214_mBq, tl208_mBq)
"""
function activity_bfd(; R_in=1.60, R_out=1.61, L=1.70)
    # Volume of cylindrical shell
    V = π * (R_out^2 - R_in^2) * L  # m³
    ρ_ptfe = 2000.0  # kg/m³
    m = V * ρ_ptfe

    a_bi214 = activity_bi214(m, "PTFE")
    a_tl208 = activity_tl208(m, "PTFE")

    return (
        mass_kg = m,
        bi214_Bq = a_bi214,
        tl208_Bq = a_tl208,
        bi214_mBq = a_bi214 * 1000,
        tl208_mBq = a_tl208 * 1000
    )
end

# ─── ICS Self-Shielding Constants ─────────────────────────────────────────────

# Copper attenuation properties at ~2.5 MeV (Bi-214: 2.448 MeV, Tl-208: 2.615 MeV)
const μ_over_ρ_Cu = 0.039  # cm²/g (mass attenuation coefficient)
const ρ_Cu = 8.96          # g/cm³
const μ_Cu = μ_over_ρ_Cu * ρ_Cu  # cm⁻¹ ≈ 0.35 cm⁻¹
const λ_Cu = 1.0 / μ_Cu    # cm ≈ 2.86 cm (attenuation length)

"""
    activity_ics(; R_in=1.61, R_out=1.64, L=1.70, L_endcap=0.03)

Compute activity of the ICS (Inner Copper Shield) - barrel + endcaps.
This is the inner 3 cm layer that contributes fully (no self-shielding).

# Arguments
- `R_in`: Inner radius (m), default 1.61 m
- `R_out`: Outer radius (m), default 1.64 m
- `L`: Barrel length (m), default 1.70 m
- `L_endcap`: Endcap thickness (m), default 0.03 m (3 cm)

# Returns
- Named tuple with barrel, endcaps, and total activities
"""
function activity_ics(; R_in=1.61, R_out=1.64, L=1.70, L_endcap=0.03)
    ρ_cu = 8960.0  # kg/m³

    # Barrel (cylindrical shell)
    V_barrel = π * (R_out^2 - R_in^2) * L
    m_barrel = V_barrel * ρ_cu

    # Endcaps (2 × solid disks)
    V_endcap = π * R_out^2 * L_endcap
    m_endcaps = 2 * V_endcap * ρ_cu

    # Total
    m_total = m_barrel + m_endcaps

    a_bi214 = activity_bi214(m_total, "Cu")
    a_tl208 = activity_tl208(m_total, "Cu")

    return (
        mass_barrel_kg = m_barrel,
        mass_endcaps_kg = m_endcaps,
        mass_total_kg = m_total,
        bi214_Bq = a_bi214,
        tl208_Bq = a_tl208,
        bi214_mBq = a_bi214 * 1000,
        tl208_mBq = a_tl208 * 1000
    )
end

"""
    activity_ics_outer_shell(; R_in=1.64, R_out=1.76, L=1.70,
                               n_shells=24, λ=λ_Cu)

Compute the effective activity contribution from the outer ICS layer (self-shielded).

The outer layer extends from 3 cm to 15 cm depth. Gammas emitted from this region
are partially absorbed before reaching the TPC. We use discrete shells to compute
the effective contribution, weighting each shell by its escape probability.

For a shell at depth d from the inner surface, the escape fraction is:
    f_escape(d) = exp(-d/λ)

# Arguments
- `R_in`: Inner radius of outer layer (m), default 1.64 m (= 1.61 + 0.03)
- `R_out`: Outer radius (m), default 1.76 m (= 1.61 + 0.15)
- `L`: Barrel length (m), default 1.70 m
- `n_shells`: Number of discrete shells for integration, default 24 (0.5 cm each)
- `λ`: Attenuation length (cm), default λ_Cu ≈ 2.86 cm

# Returns
- Named tuple with:
  - mass_total_kg: Total mass of outer layer
  - mass_effective_kg: Effective mass after self-shielding
  - shielding_factor: Fraction of activity that escapes (m_eff/m_total)
  - bi214_Bq, tl208_Bq: Effective activities
  - bi214_mBq, tl208_mBq: Activities in mBq
"""
function activity_ics_outer_shell(; R_in=1.64, R_out=1.76, L=1.70,
                                    n_shells=24, λ=λ_Cu)
    ρ_cu = 8960.0  # kg/m³

    # Shell thickness
    t_total = R_out - R_in  # Total thickness of outer layer (m)
    dr = t_total / n_shells  # Thickness of each shell (m)

    # Convert λ to meters
    λ_m = λ / 100.0  # cm → m

    m_total = 0.0
    m_effective = 0.0

    for i in 1:n_shells
        # Shell boundaries
        r_inner = R_in + (i - 1) * dr
        r_outer = R_in + i * dr
        r_mid = (r_inner + r_outer) / 2.0

        # Shell volume and mass (barrel only)
        V_shell = π * (r_outer^2 - r_inner^2) * L
        m_shell = V_shell * ρ_cu

        # Depth from inner surface of outer layer (m)
        d = r_mid - R_in

        # Escape fraction (exponential attenuation)
        f_escape = exp(-d / λ_m)

        # Accumulate
        m_total += m_shell
        m_effective += m_shell * f_escape
    end

    # Shielding factor
    shielding_factor = m_effective / m_total

    # Effective activities
    a_bi214 = activity_bi214(m_effective, "Cu")
    a_tl208 = activity_tl208(m_effective, "Cu")

    return (
        mass_total_kg = m_total,
        mass_effective_kg = m_effective,
        shielding_factor = shielding_factor,
        thickness_cm = t_total * 100,
        n_shells = n_shells,
        λ_cm = λ,
        bi214_Bq = a_bi214,
        tl208_Bq = a_tl208,
        bi214_mBq = a_bi214 * 1000,
        tl208_mBq = a_tl208 * 1000
    )
end

"""
    activity_ics_total(; R_in_inner=1.61, t_inner=0.03, t_outer=0.12,
                         L=1.70, L_endcap=0.03, n_shells=24, λ=λ_Cu)

Compute total ICS activity including both inner (unshielded) and outer (self-shielded) layers.

The ICS is 15 cm thick total:
- Inner 3 cm: contributes fully (no self-shielding)
- Outer 12 cm: partially self-shielded, contribution weighted by escape probability

# Arguments
- `R_in_inner`: Inner radius of ICS (m), default 1.61 m
- `t_inner`: Thickness of inner unshielded layer (m), default 0.03 m (3 cm)
- `t_outer`: Thickness of outer shielded layer (m), default 0.12 m (12 cm)
- `L`: Barrel length (m), default 1.70 m
- `L_endcap`: Endcap thickness (m), default 0.03 m (only inner layer for endcaps)
- `n_shells`: Number of shells for outer layer integration, default 24
- `λ`: Attenuation length (cm), default λ_Cu

# Returns
- Named tuple with inner, outer, and total contributions
"""
function activity_ics_total(; R_in_inner=1.61, t_inner=0.03, t_outer=0.12,
                              L=1.70, L_endcap=0.03, n_shells=24, λ=λ_Cu)
    # Inner layer (unshielded) - barrel + endcaps
    R_out_inner = R_in_inner + t_inner
    inner = activity_ics(; R_in=R_in_inner, R_out=R_out_inner, L=L, L_endcap=L_endcap)

    # Outer layer (self-shielded) - barrel only
    R_in_outer = R_out_inner
    R_out_outer = R_in_outer + t_outer
    outer = activity_ics_outer_shell(; R_in=R_in_outer, R_out=R_out_outer, L=L,
                                       n_shells=n_shells, λ=λ)

    # Total
    m_total = inner.mass_total_kg + outer.mass_total_kg
    m_effective = inner.mass_total_kg + outer.mass_effective_kg

    bi214_total = inner.bi214_Bq + outer.bi214_Bq
    tl208_total = inner.tl208_Bq + outer.tl208_Bq

    return (
        # Inner layer (3 cm, unshielded)
        inner_mass_kg = inner.mass_total_kg,
        inner_bi214_mBq = inner.bi214_mBq,
        inner_tl208_mBq = inner.tl208_mBq,
        # Outer layer (12 cm, self-shielded)
        outer_mass_total_kg = outer.mass_total_kg,
        outer_mass_effective_kg = outer.mass_effective_kg,
        outer_shielding_factor = outer.shielding_factor,
        outer_bi214_mBq = outer.bi214_mBq,
        outer_tl208_mBq = outer.tl208_mBq,
        # Total ICS
        mass_total_kg = m_total,
        mass_effective_kg = m_effective,
        bi214_Bq = bi214_total,
        tl208_Bq = tl208_total,
        bi214_mBq = bi214_total * 1000,
        tl208_mBq = tl208_total * 1000
    )
end

"""
    print_ics_shielding_summary(; R_in_inner=1.61, t_inner=0.03, t_outer=0.12,
                                  L=1.70, L_endcap=0.03, n_shells=24, λ=λ_Cu)

Print a detailed summary of ICS activity with self-shielding calculations.
"""
function print_ics_shielding_summary(; R_in_inner=1.61, t_inner=0.03, t_outer=0.12,
                                       L=1.70, L_endcap=0.03, n_shells=24, λ=λ_Cu)
    result = activity_ics_total(; R_in_inner, t_inner, t_outer, L, L_endcap, n_shells, λ)

    println("\n" * "="^70)
    println("          ICS Self-Shielding Activity Summary")
    println("="^70)

    println("\n── Geometry ──")
    @printf("  Inner radius:          R_in = %.2f m\n", R_in_inner)
    @printf("  Inner layer thickness: t_inner = %.0f cm (unshielded)\n", t_inner * 100)
    @printf("  Outer layer thickness: t_outer = %.0f cm (self-shielded)\n", t_outer * 100)
    @printf("  Total thickness:       t_total = %.0f cm\n", (t_inner + t_outer) * 100)
    @printf("  Barrel length:         L = %.2f m\n", L)

    println("\n── Attenuation Properties (Cu at ~2.5 MeV) ──")
    @printf("  μ/ρ = %.3f cm²/g\n", μ_over_ρ_Cu)
    @printf("  ρ   = %.2f g/cm³\n", ρ_Cu)
    @printf("  μ   = %.3f cm⁻¹\n", μ_Cu)
    @printf("  λ   = %.2f cm (attenuation length)\n", λ)

    println("\n── Inner Layer (3 cm, unshielded) ──")
    @printf("  Mass:       %.2f kg\n", result.inner_mass_kg)
    @printf("  Bi-214:     %.4f mBq\n", result.inner_bi214_mBq)
    @printf("  Tl-208:     %.4f mBq\n", result.inner_tl208_mBq)

    println("\n── Outer Layer (12 cm, self-shielded) ──")
    @printf("  Total mass:     %.2f kg\n", result.outer_mass_total_kg)
    @printf("  Effective mass: %.2f kg (after shielding)\n", result.outer_mass_effective_kg)
    @printf("  Shielding factor: %.3f (%.1f%% escapes)\n",
            result.outer_shielding_factor, result.outer_shielding_factor * 100)
    @printf("  Bi-214:     %.4f mBq (effective)\n", result.outer_bi214_mBq)
    @printf("  Tl-208:     %.4f mBq (effective)\n", result.outer_tl208_mBq)

    println("\n── Total ICS (inner + outer) ──")
    @printf("  Total mass:     %.2f kg\n", result.mass_total_kg)
    @printf("  Effective mass: %.2f kg\n", result.mass_effective_kg)
    @printf("  Bi-214:     %.4f mBq\n", result.bi214_mBq)
    @printf("  Tl-208:     %.4f mBq\n", result.tl208_mBq)

    # Compare with inner-only
    inner_only = activity_ics(; R_in=R_in_inner, R_out=R_in_inner+t_inner, L=L, L_endcap=L_endcap)
    println("\n── Comparison ──")
    @printf("  Inner 3cm only:  Bi-214 = %.4f mBq, Tl-208 = %.4f mBq\n",
            inner_only.bi214_mBq, inner_only.tl208_mBq)
    @printf("  With outer 12cm: Bi-214 = %.4f mBq, Tl-208 = %.4f mBq\n",
            result.bi214_mBq, result.tl208_mBq)
    @printf("  Increase factor: Bi-214 = %.2fx, Tl-208 = %.2fx\n",
            result.bi214_mBq / inner_only.bi214_mBq,
            result.tl208_mBq / inner_only.tl208_mBq)

    println("\n" * "="^70)

    return result
end

"""
    activity_sipm_boards(; D=3.20, t=0.2e-3)

Compute activity of the SiPM Kapton boards in the DSP.

The SiPM boards are Kapton disks holding the SiPM arrays.

# Arguments
- `D`: Diameter of the board (m), default 3.20 m
- `t`: Thickness of the Kapton board (m), default 0.2 mm

# Returns
- Named tuple with mass and activities
"""
function activity_sipm_boards(; D=3.20, t=0.2e-3)
    ρ_kapton = 1420.0  # kg/m³

    # Disk volume and mass
    R = D / 2.0
    V = π * R^2 * t  # m³
    m = V * ρ_kapton  # kg

    a_bi214 = activity_bi214(m, "Kapton")
    a_tl208 = activity_tl208(m, "Kapton")

    return (
        diameter_m = D,
        thickness_mm = t * 1000,
        mass_kg = m,
        bi214_Bq = a_bi214,
        tl208_Bq = a_tl208,
        bi214_mBq = a_bi214 * 1000,
        tl208_mBq = a_tl208 * 1000
    )
end

"""
    activity_dsp(; D=2.30, L=0.10, d_channel=6e-3, pitch=10e-3)

Compute activity of the DSP (Dense Silicon Plane) - HDPE honeycomb.

Uses the honeycomb mass calculation function.

# Returns
- Named tuple with mass and activities
"""
function activity_dsp(; D=2.30, L=0.10, d_channel=6e-3, pitch=10e-3)
    result = mass_dsp_honeycomb(; D, L, d_channel, pitch, ρ=ρ_HDPE)
    m = result.mass_honeycomb

    a_bi214 = activity_bi214(m, "HDPE")
    a_tl208 = activity_tl208(m, "HDPE")

    return (
        mass_kg = m,
        fill_factor = result.fill_factor,
        bi214_Bq = a_bi214,
        tl208_Bq = a_tl208,
        bi214_mBq = a_bi214 * 1000,
        tl208_mBq = a_tl208 * 1000
    )
end

"""
    activity_cathode(; D=3.20, w=w_mesh_default, M=M_mesh_default)

Compute activity of the cathode mesh (Fe316Ti).

Uses the cathode mesh mass calculation function.

# Returns
- Named tuple with mass and activities
"""
function activity_cathode(; D=3.20, w=w_mesh_default, M=M_mesh_default)
    result = mass_cathode_mesh(; D, w, M, ρ=ρ_Fe316Ti)
    m = result.mass

    a_bi214 = activity_bi214(m, "Fe316Ti")
    a_tl208 = activity_tl208(m, "Fe316Ti")

    return (
        mass_kg = m,
        bi214_Bq = a_bi214,
        tl208_Bq = a_tl208,
        bi214_mBq = a_bi214 * 1000,
        tl208_mBq = a_tl208 * 1000
    )
end

"""
    activity_mars(; N=3, d_o=5e-3, t_w=0.5e-3, R=1.6,
                    perimeter=0.165, t_skin=50e-6,
                    belt_w=10e-3, belt_t=125e-6,
                    cable_w=10e-3, cable_t=50e-6,
                    n_ribs=10, m_rib=5e-3,
                    Lx=0.16, Ly=0.16, Lz=0.02, t_p=0.002, m_extra=0.03)

Compute activity of the MARS system (composite: Ti longerons + Kapton + HDPE).

Components:
- Ti Gr.5 longerons
- Kapton skin, belt, cables
- HDPE ribs and ion carrier plates

# Returns
- Named tuple with per-component and total activities
"""
function activity_mars(; N=3, d_o=5e-3, t_w=0.5e-3, R=1.6,
                         perimeter=0.165, t_skin=50e-6,
                         belt_w=10e-3, belt_t=125e-6,
                         cable_w=10e-3, cable_t=50e-6,
                         n_ribs=10, m_rib=5e-3,
                         Lx=0.16, Ly=0.16, Lz=0.02, t_p=0.002, m_extra=0.0)
    # Masses (note: m_extra excluded from activity calc as it's CMOS/screws)
    m_long = mass_longerons(; N, d_o, t_w, R)
    m_skin = mass_mars_skin(; perimeter, t_skin, R)
    m_belt = mass_mars_belt(; width=belt_w, thickness=belt_t, R)
    m_cables = mass_mars_cables(; width=cable_w, thickness=cable_t, R)
    m_ribs = mass_mars_ribs(; n_ribs, m_rib)
    m_plates = mass_ion_carrier(; Lx, Ly, Lz, t_p, m_extra=0.0)  # HDPE only

    # Activity by material group
    # Ti: longerons
    a_bi214_ti = activity_bi214(m_long, "Ti")
    a_tl208_ti = activity_tl208(m_long, "Ti")

    # Kapton: skin + belt + cables
    m_kapton = m_skin + m_belt + m_cables
    a_bi214_kapton = activity_bi214(m_kapton, "Kapton")
    a_tl208_kapton = activity_tl208(m_kapton, "Kapton")

    # HDPE: ribs + plates
    m_hdpe = m_ribs + m_plates
    a_bi214_hdpe = activity_bi214(m_hdpe, "HDPE")
    a_tl208_hdpe = activity_tl208(m_hdpe, "HDPE")

    # Totals
    m_total = m_long + m_kapton + m_hdpe
    a_bi214_total = a_bi214_ti + a_bi214_kapton + a_bi214_hdpe
    a_tl208_total = a_tl208_ti + a_tl208_kapton + a_tl208_hdpe

    return (
        # Masses
        mass_ti_kg = m_long,
        mass_kapton_kg = m_kapton,
        mass_hdpe_kg = m_hdpe,
        mass_total_kg = m_total,
        # Per-material activities (Bq)
        bi214_ti_Bq = a_bi214_ti,
        bi214_kapton_Bq = a_bi214_kapton,
        bi214_hdpe_Bq = a_bi214_hdpe,
        tl208_ti_Bq = a_tl208_ti,
        tl208_kapton_Bq = a_tl208_kapton,
        tl208_hdpe_Bq = a_tl208_hdpe,
        # Total activities
        bi214_Bq = a_bi214_total,
        tl208_Bq = a_tl208_total,
        bi214_mBq = a_bi214_total * 1000,
        tl208_mBq = a_tl208_total * 1000
    )
end

# ─── Total Activity Summary ────────────────────────────────────────────────

"""
    activity_itaca_summary()

Compute and return total activity summary for all ITACA detector components.

# Returns
- Named tuple with all component activities and totals
"""
function activity_itaca_summary()
    bfd = activity_bfd()
    ics = activity_ics()
    dsp = activity_dsp()
    sipm_boards = activity_sipm_boards()
    cathode = activity_cathode()
    mars = activity_mars()

    total_bi214 = bfd.bi214_Bq + ics.bi214_Bq + dsp.bi214_Bq + sipm_boards.bi214_Bq + cathode.bi214_Bq + mars.bi214_Bq
    total_tl208 = bfd.tl208_Bq + ics.tl208_Bq + dsp.tl208_Bq + sipm_boards.tl208_Bq + cathode.tl208_Bq + mars.tl208_Bq
    total_mass = bfd.mass_kg + ics.mass_total_kg + dsp.mass_kg + sipm_boards.mass_kg + cathode.mass_kg + mars.mass_total_kg

    return (
        bfd = bfd,
        ics = ics,
        dsp = dsp,
        sipm_boards = sipm_boards,
        cathode = cathode,
        mars = mars,
        total_mass_kg = total_mass,
        total_bi214_Bq = total_bi214,
        total_tl208_Bq = total_tl208,
        total_bi214_mBq = total_bi214 * 1000,
        total_tl208_mBq = total_tl208 * 1000
    )
end

"""
    print_activity_summary()

Print a formatted summary of all ITACA detector component activities.
"""
function print_activity_summary()
    s = activity_itaca_summary()

    println("\n" * "="^70)
    println("          ITACA Detector Activity Summary")
    println("="^70)

    println("\n── Component Activities ──\n")
    println("  Component          Material     Mass (kg)    Bi-214 (mBq)   Tl-208 (mBq)")
    println("  " * "-"^66)

    @printf("  BFD                PTFE         %8.2f     %10.4f     %10.4f\n",
            s.bfd.mass_kg, s.bfd.bi214_mBq, s.bfd.tl208_mBq)

    @printf("  ICS (barrel+caps)  Cu           %8.2f     %10.4f     %10.4f\n",
            s.ics.mass_total_kg, s.ics.bi214_mBq, s.ics.tl208_mBq)

    @printf("  DSP (honeycomb)    HDPE         %8.2f     %10.4f     %10.4f\n",
            s.dsp.mass_kg, s.dsp.bi214_mBq, s.dsp.tl208_mBq)

    @printf("  SiPM boards        Kapton       %8.4f     %10.4f     %10.4f\n",
            s.sipm_boards.mass_kg, s.sipm_boards.bi214_mBq, s.sipm_boards.tl208_mBq)

    @printf("  Cathode (mesh)     Fe316Ti      %8.4f     %10.6f     %10.6f\n",
            s.cathode.mass_kg, s.cathode.bi214_mBq, s.cathode.tl208_mBq)

    println("\n── MARS System (composite) ──\n")
    @printf("    Ti longerons                  %8.4f     %10.6f     %10.6f\n",
            s.mars.mass_ti_kg, s.mars.bi214_ti_Bq * 1000, s.mars.tl208_ti_Bq * 1000)
    @printf("    Kapton (skin+belt+cables)     %8.4f     %10.6f     %10.6f\n",
            s.mars.mass_kapton_kg, s.mars.bi214_kapton_Bq * 1000, s.mars.tl208_kapton_Bq * 1000)
    @printf("    HDPE (ribs+plates)            %8.4f     %10.6f     %10.6f\n",
            s.mars.mass_hdpe_kg, s.mars.bi214_hdpe_Bq * 1000, s.mars.tl208_hdpe_Bq * 1000)
    @printf("  MARS Total                      %8.4f     %10.6f     %10.6f\n",
            s.mars.mass_total_kg, s.mars.bi214_mBq, s.mars.tl208_mBq)

    println("\n  " * "-"^66)
    @printf("  TOTAL              ---         %8.2f     %10.4f     %10.4f\n",
            s.total_mass_kg, s.total_bi214_mBq, s.total_tl208_mBq)

    println("\n── Specific Activities Used (Bq/kg) ──\n")
    println("  Material      Bi-214         Tl-208")
    println("  " * "-"^40)
    @printf("  Fe316Ti       %.2e     %.2e\n", SA_Bi214["Fe316Ti"], SA_Tl208["Fe316Ti"])
    @printf("  Cu            %.2e     %.2e\n", SA_Bi214["Cu"], SA_Tl208["Cu"])
    @printf("  PTFE          %.2e     %.2e\n", SA_Bi214["PTFE"], SA_Tl208["PTFE"])
    @printf("  HDPE          %.2e     %.2e\n", SA_Bi214["HDPE"], SA_Tl208["HDPE"])
    @printf("  Kapton        %.2e     %.2e\n", SA_Bi214["Kapton"], SA_Tl208["Kapton"])
    @printf("  Ti            %.2e     %.2e\n", SA_Bi214["Ti"], SA_Tl208["Ti"])

    println("\n" * "="^70)

    return s
end

# ─── DataFrame Output Functions (Standalone) ─────────────────────────────────

"""
    create_activity_dataframe()

Create a DataFrame with activity summary for all ITACA detector components.
Uses the standalone activity functions (no JHD5t detector object required).

Includes two ICS entries:
- ICS_3cm: Inner 3 cm layer (unshielded, full contribution)
- ICS_15cm: Full 15 cm layer (inner 3 cm + outer 12 cm with self-shielding)

# Returns
- DataFrame with columns: Component, Material, Mass_kg, Bi214_mBq, Tl208_mBq
"""
function create_activity_dataframe()
    s = activity_itaca_summary()

    # Get ICS with self-shielding for full 15 cm
    ics_full = activity_ics_total()

    # Build arrays for DataFrame columns
    components = String[]
    materials = String[]
    masses = Float64[]
    bi214_mBq = Float64[]
    tl208_mBq = Float64[]

    # BFD
    push!(components, "BFD")
    push!(materials, "PTFE")
    push!(masses, s.bfd.mass_kg)
    push!(bi214_mBq, s.bfd.bi214_mBq)
    push!(tl208_mBq, s.bfd.tl208_mBq)

    # ICS - Inner 3 cm only (unshielded)
    push!(components, "ICS_3cm")
    push!(materials, "Cu")
    push!(masses, s.ics.mass_total_kg)
    push!(bi214_mBq, s.ics.bi214_mBq)
    push!(tl208_mBq, s.ics.tl208_mBq)

    # ICS - Full 15 cm with self-shielding (inner 3cm + outer 12cm shielded)
    push!(components, "ICS_15cm")
    push!(materials, "Cu")
    push!(masses, ics_full.mass_effective_kg)  # Effective mass after shielding
    push!(bi214_mBq, ics_full.bi214_mBq)
    push!(tl208_mBq, ics_full.tl208_mBq)

    # DSP
    push!(components, "DSP")
    push!(materials, "HDPE")
    push!(masses, s.dsp.mass_kg)
    push!(bi214_mBq, s.dsp.bi214_mBq)
    push!(tl208_mBq, s.dsp.tl208_mBq)

    # SiPM boards (Kapton)
    push!(components, "SiPM_boards")
    push!(materials, "Kapton")
    push!(masses, s.sipm_boards.mass_kg)
    push!(bi214_mBq, s.sipm_boards.bi214_mBq)
    push!(tl208_mBq, s.sipm_boards.tl208_mBq)

    # Cathode
    push!(components, "Cathode")
    push!(materials, "Fe316Ti")
    push!(masses, s.cathode.mass_kg)
    push!(bi214_mBq, s.cathode.bi214_mBq)
    push!(tl208_mBq, s.cathode.tl208_mBq)

    # MARS - Ti longerons
    push!(components, "MARS_Ti")
    push!(materials, "Ti")
    push!(masses, s.mars.mass_ti_kg)
    push!(bi214_mBq, s.mars.bi214_ti_Bq * 1000)
    push!(tl208_mBq, s.mars.tl208_ti_Bq * 1000)

    # MARS - Kapton
    push!(components, "MARS_Kapton")
    push!(materials, "Kapton")
    push!(masses, s.mars.mass_kapton_kg)
    push!(bi214_mBq, s.mars.bi214_kapton_Bq * 1000)
    push!(tl208_mBq, s.mars.tl208_kapton_Bq * 1000)

    # MARS - HDPE
    push!(components, "MARS_HDPE")
    push!(materials, "HDPE")
    push!(masses, s.mars.mass_hdpe_kg)
    push!(bi214_mBq, s.mars.bi214_hdpe_Bq * 1000)
    push!(tl208_mBq, s.mars.tl208_hdpe_Bq * 1000)

    # Create DataFrame
    df = DataFrame(
        Component = components,
        Material = materials,
        Mass_kg = masses,
        Bi214_mBq = bi214_mBq,
        Tl208_mBq = tl208_mBq
    )

    # Calculate totals for both ICS scenarios
    # Total with ICS_3cm (original calculation)
    total_mass_3cm = s.total_mass_kg
    total_bi214_3cm = s.total_bi214_mBq
    total_tl208_3cm = s.total_tl208_mBq

    # Total with ICS_15cm (replace ICS_3cm contribution with ICS_15cm)
    total_mass_15cm = s.total_mass_kg - s.ics.mass_total_kg + ics_full.mass_effective_kg
    total_bi214_15cm = s.total_bi214_mBq - s.ics.bi214_mBq + ics_full.bi214_mBq
    total_tl208_15cm = s.total_tl208_mBq - s.ics.tl208_mBq + ics_full.tl208_mBq

    # Add totals rows
    totals_rows = DataFrame(
        Component = ["TOTAL_ICS3cm", "TOTAL_ICS15cm"],
        Material = ["---", "---"],
        Mass_kg = [total_mass_3cm, total_mass_15cm],
        Bi214_mBq = [total_bi214_3cm, total_bi214_15cm],
        Tl208_mBq = [total_tl208_3cm, total_tl208_15cm]
    )

    df = vcat(df, totals_rows)

    return df
end

"""
    save_activity_dataframe(; filename="itaca_activity_summary.csv", output_dir=".")

Create and save the activity summary DataFrame to a CSV file.

# Arguments
- `filename`: Name of the CSV file (default: "itaca_activity_summary.csv")
- `output_dir`: Output directory (default: current directory)

# Returns
- The DataFrame that was saved
"""
function save_activity_dataframe(; filename="itaca_activity_summary.csv", output_dir=".")
    df = create_activity_dataframe()

    # Ensure output directory exists
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Save to CSV
    filepath = joinpath(output_dir, filename)
    CSV.write(filepath, df)

    println("Activity summary saved to: $filepath")
    return df
end

# ─── Plotting Functions ───────────────────────────────────────────────────────

"""
    plot_itaca_projections(; save_plots=false, output_dir=".")

Draw XY and XZ projections of the itaca detector showing different materials with distinct colors.
Uses geometry from module constants (standalone, no detector object required).
"""
function plot_itaca_projections(; save_plots=false, output_dir=".")

    # Geometry parameters from constants
    R_fid = D_fid_default * 100.0 / 2.0  # 160 cm
    L_tpc = L_fid_default * 100.0        # 150 cm
    L_dsp = DSP_LG_FATGEM_cm             # 10 cm
    L_mrs = MRS_thickness_cm             # 10 cm
    L_internal = L_tpc + L_dsp + L_mrs   # 170 cm

    R_bfd_in = R_fid                     # 160 cm
    R_bfd_out = R_fid + BFD_thickness_cm # 161 cm
    R_ics_in = R_bfd_out                 # 161 cm
    R_ics_out = R_ics_in + ICS_thickness_cm  # 164 cm

    z_tpc_top = L_tpc / 2.0              # +75 cm
    z_tpc_bot = -L_tpc / 2.0             # -75 cm
    z_dsp_top = z_tpc_top + L_dsp        # +85 cm
    z_mrs_bot = z_tpc_bot - L_mrs        # -85 cm
    z_endcap_top = L_internal / 2.0 + ICS_thickness_cm  # +88 cm
    z_endcap_bot = -z_endcap_top         # -88 cm

    # Define material colors
    colors = Dict(
        "copper" => :orange,     # Copper shielding
        "ptfe" => :purple,       # PTFE / BFD
        "steel" => :gray,        # Central cathode
        "xenon" => :lightblue,   # Gas volume
        "dsp" => :green          # DSP region
    )

    # Create XY projection (top view)
    xy_plot = plot(
        title="itaca Detector - XY Projection (Top View)",
        xlabel="X (cm)", ylabel="Y (cm)",
        aspect_ratio=:equal,
        legend=:topright,
        size=(800, 800)
    )

    # Plot concentric circles
    θ = range(0, 2π, length=100)

    # ICS (Copper shield outer)
    x_ics_out = R_ics_out .* cos.(θ)
    y_ics_out = R_ics_out .* sin.(θ)
    x_ics_in = R_ics_in .* cos.(θ)
    y_ics_in = R_ics_in .* sin.(θ)

    plot!(xy_plot, x_ics_out, y_ics_out, linewidth=3, color=colors["copper"],
          label="ICS (Copper)", fill=false)
    plot!(xy_plot, x_ics_in, y_ics_in, linewidth=2, color=colors["copper"],
          label="", fill=false, alpha=0.7)

    # BFD (PTFE)
    x_bfd_out = R_bfd_out .* cos.(θ)
    y_bfd_out = R_bfd_out .* sin.(θ)
    x_bfd_in = R_bfd_in .* cos.(θ)
    y_bfd_in = R_bfd_in .* sin.(θ)

    plot!(xy_plot, x_bfd_out, y_bfd_out, linewidth=2, color=colors["ptfe"],
          label="BFD (PTFE)", fill=false)
    plot!(xy_plot, x_bfd_in, y_bfd_in, linewidth=1, color=colors["ptfe"],
          label="", fill=false, alpha=0.7)

    # Gas volume (TPC)
    x_gas = R_fid .* cos.(θ)
    y_gas = R_fid .* sin.(θ)
    plot!(xy_plot, x_gas, y_gas, linewidth=1, color=colors["xenon"],
          label="TPC (Xenon)", fill=true, alpha=0.3)

    # Central cathode (visible as a cross in XY view)
    plot!(xy_plot, [-R_fid, R_fid], [0, 0], linewidth=2, color=colors["steel"],
          label="Cathode")
    plot!(xy_plot, [0, 0], [-R_fid, R_fid], linewidth=2, color=colors["steel"],
          label="")

    # Add center point
    scatter!(xy_plot, [0], [0], color=:red, markersize=3, label="Center")

    # Create XZ projection (side view)
    xz_plot = plot(
        title="itaca Detector - XZ Projection (Side View)",
        xlabel="X (cm)", ylabel="Z (cm)",
        aspect_ratio=:equal,
        legend=:topright,
        size=(1000, 600)
    )

    # ICS barrel (copper)
    plot!(xz_plot, [R_ics_out, R_ics_out], [-L_internal/2, L_internal/2], linewidth=3,
          color=colors["copper"], label="ICS (Copper)")
    plot!(xz_plot, [-R_ics_out, -R_ics_out], [-L_internal/2, L_internal/2], linewidth=3,
          color=colors["copper"], label="")
    plot!(xz_plot, [R_ics_in, R_ics_in], [-L_internal/2, L_internal/2], linewidth=2,
          color=colors["copper"], label="", alpha=0.7)
    plot!(xz_plot, [-R_ics_in, -R_ics_in], [-L_internal/2, L_internal/2], linewidth=2,
          color=colors["copper"], label="", alpha=0.7)

    # ICS endcaps (copper)
    plot!(xz_plot, [-R_ics_out, R_ics_out], [z_endcap_top, z_endcap_top], linewidth=3,
          color=colors["copper"], label="")
    plot!(xz_plot, [-R_ics_out, R_ics_out], [z_endcap_bot, z_endcap_bot], linewidth=3,
          color=colors["copper"], label="")
    plot!(xz_plot, [-R_ics_out, R_ics_out], [L_internal/2, L_internal/2], linewidth=2,
          color=colors["copper"], label="", alpha=0.7)
    plot!(xz_plot, [-R_ics_out, R_ics_out], [-L_internal/2, -L_internal/2], linewidth=2,
          color=colors["copper"], label="", alpha=0.7)

    # BFD barrel (PTFE)
    plot!(xz_plot, [R_bfd_out, R_bfd_out], [-L_internal/2, L_internal/2], linewidth=2,
          color=colors["ptfe"], label="BFD (PTFE)")
    plot!(xz_plot, [-R_bfd_out, -R_bfd_out], [-L_internal/2, L_internal/2], linewidth=2,
          color=colors["ptfe"], label="")
    plot!(xz_plot, [R_bfd_in, R_bfd_in], [-L_internal/2, L_internal/2], linewidth=1,
          color=colors["ptfe"], label="", alpha=0.7)
    plot!(xz_plot, [-R_bfd_in, -R_bfd_in], [-L_internal/2, L_internal/2], linewidth=1,
          color=colors["ptfe"], label="", alpha=0.7)

    # DSP region (HDPE honeycomb, anode side)
    plot!(xz_plot, [-R_fid, R_fid, R_fid, -R_fid, -R_fid],
          [z_tpc_top, z_tpc_top, z_dsp_top, z_dsp_top, z_tpc_top],
          linewidth=1, color=colors["dsp"], label="DSP",
          fill=true, alpha=0.3)

    # MRS region (MARS region, cathode side)
    plot!(xz_plot, [-R_fid, R_fid, R_fid, -R_fid, -R_fid],
          [z_mrs_bot, z_mrs_bot, z_tpc_bot, z_tpc_bot, z_mrs_bot],
          linewidth=1, color=colors["xenon"], label="MRS",
          fill=true, alpha=0.2)

    # TPC gas volume (filled)
    plot!(xz_plot, [-R_fid, R_fid, R_fid, -R_fid, -R_fid],
          [z_tpc_bot, z_tpc_bot, z_tpc_top, z_tpc_top, z_tpc_bot],
          linewidth=1, color=colors["xenon"], label="TPC (Xenon)",
          fill=true, alpha=0.3)

    # Central cathode
    plot!(xz_plot, [-R_fid, R_fid], [z_tpc_bot, z_tpc_bot], linewidth=3,
          color=colors["steel"], label="Cathode")

    # Add center line
    plot!(xz_plot, [0, 0], [z_endcap_bot-10, z_endcap_top+10],
          color=:red, linestyle=:dash, alpha=0.5, label="Detector Axis")

    # Save plots if requested
    if save_plots
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        savefig(xy_plot, joinpath(output_dir, "itaca_xy_projection.png"))
        savefig(xz_plot, joinpath(output_dir, "itaca_xz_projection.png"))
        println("Plots saved to: $output_dir")
    end

    return (xy_plot, xz_plot)
end

"""
    plot_itaca_combined(; save_plot=false, output_dir=".")

Create a combined view of the itaca detector projections.
"""
function plot_itaca_combined(; save_plot=false, output_dir=".")
    xy_plot, xz_plot = plot_itaca_projections(save_plots=false)

    combined_plot = plot(xy_plot, xz_plot,
                        layout=(1,2),
                        size=(1600, 800),
                        plot_title="itaca Detector - Multiple Projections")

    if save_plot
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        savefig(combined_plot, joinpath(output_dir, "itaca_combined.png"))
        println("Combined plot saved to: $output_dir")
    end

    return combined_plot
end

# ─── Verbose Print Functions ──────────────────────────────────────────────────

"""
    print_verbose_summary()

Print verbose summary including all component details: DSP honeycomb, cathode mesh,
and MARS system breakdown.
"""
function print_verbose_summary()
    # Print DSP honeycomb summary
    print_dsp_honeycomb_summary()

    # Print cathode mesh summary
    print_cathode_mesh_summary()

    # Print MARS summary (from mars_pars.jl)
    mars_full_summary()

    # Print ICS self-shielding summary
    print_ics_shielding_summary()

    # Print activity summary
    print_activity_summary()
end

# ─── Command Line Argument Parsing ────────────────────────────────────────────

function parse_commandline()
    s = ArgParseSettings(
        description = "ITACA Detector Activity Calculator",
        version = "1.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--plots"
            help = "Create and save detector projection plots (default: false)"
            action = :store_true
        "--no-dataframe"
            help = "Disable saving activity summary DataFrame (default: save)"
            action = :store_true
        "--print"
            help = "Print mode: 'concise' (activity summary only) or 'verbose' (all details)"
            arg_type = String
            default = "concise"
        "--output-dir"
            help = "Output directory for files"
            arg_type = String
            default = "."
    end

    return parse_args(s)
end

# ─── Main Function ────────────────────────────────────────────────────────────

function main()
    # Parse command line arguments
    args = parse_commandline()

    # Determine effective dataframe flag (--no-dataframe overrides --dataframe)
    save_df = !args["no-dataframe"]

    # Get print mode
    print_mode = lowercase(args["print"])
    output_dir = args["output-dir"]

    println("ITACA Detector Analysis")
    println("="^60)

    # Print based on mode
    if print_mode == "verbose"
        print_verbose_summary()
    else
        # Concise mode: only activity summary
        print_activity_summary()
    end

    # Save DataFrame if requested
    if save_df
        println("\nSaving activity DataFrame...")
        df = save_activity_dataframe(output_dir=output_dir)
    end

    # Create plots if requested
    if args["plots"]
        println("\nCreating detector visualizations...")
        plots_dir = joinpath(output_dir, "itaca_plots")
        xy_plot, xz_plot = plot_itaca_projections(save_plots=true, output_dir=plots_dir)
        combined_plot = plot_itaca_combined(save_plot=true, output_dir=plots_dir)

        println("\nPlots saved to: $plots_dir/")
        println("   - itaca_xy_projection.png")
        println("   - itaca_xz_projection.png")
        println("   - itaca_combined.png")
    end

    println("\nITACA detector analysis completed!")
end

# Run the script if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end