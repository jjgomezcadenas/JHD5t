#!/usr/bin/env julia

"""
ITACA Detector Functions
========================
Core functions for ITACA detector calculations including:
- Fiducial mass and geometry
- Ion and electron diffusion
- Pressure vessel sizing (ASME)
- MARS kinematics
- Electrode voltages
- Structural analysis
"""

using Printf
using Plots

# =============================================================================
# ITACA Default Parameters
# =============================================================================

# Fiducial dimensions
const L_fid_default = 1.5   # m (fiducial length)
const D_fid_default = 3.2   # m (fiducial diameter)

# Electric field and drift
const E_default = 300.0     # V/cm (drift field)
const v_drift_default = 15.0  # cm/s (ion drift velocity at E=300 V/cm, P=15 bar)
const v_electron_mm_per_us = 1.0  # mm/μs (electron drift velocity)

# Voltages
const V_cathode = 0.5       # kV (cathode voltage)
const V_EL = 15.0           # kV (EL/FAT-GEM voltage)

# Xenon density at 15 bar (kg/m³)
const ρ_Xe_15bar = 87.0

# Physical constants
const k_B = 1.380649e-23   # Boltzmann constant (J/K)
const q_e = 1.602176634e-19 # Elementary charge (C)

# FAT-GEM voltage parameters (same as V_EL)
const V_FATGEM = V_EL  # kV across FAT-GEM

# Cathode mesh default parameters
const w_mesh_default = 200e-6  # m (wire diameter, 200 μm)
const M_mesh_default = 5e-3    # m (mesh pitch, 5 mm)

# =============================================================================
# Cathode Mesh Functions
# =============================================================================

"""
    mesh_transparency(w::Real, M::Real)

Compute the optical/geometric transparency of a square wire mesh.

For a square mesh with wire diameter w and pitch M, the transparency is:
    f = (1 - w/M)²

# Arguments
- `w`: Wire diameter in meters
- `M`: Mesh pitch (wire spacing) in meters

# Returns
- Transparency f (dimensionless, 0 to 1)
"""
function mesh_transparency(w::Real, M::Real)
    return (1 - w/M)^2
end

"""
    mesh_transparency_um_mm(w_um::Real, M_mm::Real)

Same as `mesh_transparency` but with w in μm and M in mm.

# Arguments
- `w_um`: Wire diameter in μm
- `M_mm`: Mesh pitch in mm

# Returns
- Transparency f (dimensionless, 0 to 1)
"""
function mesh_transparency_um_mm(w_um::Real, M_mm::Real)
    w_m = w_um * 1e-6
    M_m = M_mm * 1e-3
    return mesh_transparency(w_m, M_m)
end

# =============================================================================
# Electron and Ion Drift Time Functions
# =============================================================================

"""
    electron_drift_time_ms(L_cm::Real; v_e::Real=v_electron_mm_per_us)

Compute the electron drift time for a given drift length.

# Arguments
- `L_cm`: Drift length in cm
- `v_e`: Electron drift velocity in mm/μs (default: 1.0 mm/μs)

# Returns
- Drift time in milliseconds
"""
function electron_drift_time_ms(L_cm::Real; v_e::Real=v_electron_mm_per_us)
    L_mm = L_cm * 10.0  # cm → mm
    # v_e is in mm/μs, so t = L/v gives time in μs
    t_us = L_mm / v_e
    return t_us / 1000.0  # μs → ms
end

"""
    ion_drift_time_s(L_cm::Real; v_ion::Real=v_drift_default)

Compute the ion drift time for a given drift length.

# Arguments
- `L_cm`: Drift length in cm
- `v_ion`: Ion drift velocity in cm/s (default: 15 cm/s)

# Returns
- Drift time in seconds
"""
function ion_drift_time_s(L_cm::Real; v_ion::Real=v_drift_default)
    return L_cm / v_ion
end

"""
    ion_drift_velocity_at_field(P_bar::Real, E::Real; μ_stp::Real=0.74)

Compute the ion drift velocity at a given pressure and electric field.

Uses the relation: v = μ(P) × E, where μ(P) = μ_STP × (P_STP / P)

# Arguments
- `P_bar`: Pressure in bar
- `E`: Electric field in V/cm
- `μ_stp`: Reduced mobility at STP in cm²/V/s (default: 0.74 for Xe₂⁺)

# Returns
- Drift velocity in cm/s
"""
function ion_drift_velocity_at_field(P_bar::Real, E::Real; μ_stp::Real=0.74)
    μ_P = μ_stp / P_bar  # Mobility scales as 1/P
    return μ_P * E
end

# =============================================================================
# ASME Pressure Vessel Functions (Section VIII, Division 1)
# =============================================================================

"""
    cylinder_wall_thickness(P::Real, R_i::Real, S::Real, E_w::Real)

Compute the minimum wall thickness for a thin-walled cylindrical shell
using the ASME formula (Section VIII, Division 1):

    t = (P × R_i) / (S × E_w - 0.6 × P)

# Arguments
- `P`: Internal design pressure (MPa)
- `R_i`: Inner radius (mm)
- `S`: Maximum allowable stress of the material (MPa)
- `E_w`: Joint efficiency factor (dimensionless, 0 < E_w ≤ 1)

# Returns
- Wall thickness t in mm
"""
function cylinder_wall_thickness(P::Real, R_i::Real, S::Real, E_w::Real)
    return (P * R_i) / (S * E_w - 0.6 * P)
end

"""
    hemisphere_wall_thickness(P::Real, R_i::Real, S::Real, E_w::Real)

Compute the minimum wall thickness for a hemispherical head
using the ASME formula (Section VIII, Division 1):

    t_h = (P × R_i) / (2 × S × E_w - 0.2 × P)

# Arguments
- `P`: Internal design pressure (MPa)
- `R_i`: Inner radius (mm)
- `S`: Maximum allowable stress of the material (MPa)
- `E_w`: Joint efficiency factor (dimensionless, 0 < E_w ≤ 1)

# Returns
- Wall thickness t_h in mm
"""
function hemisphere_wall_thickness(P::Real, R_i::Real, S::Real, E_w::Real)
    return (P * R_i) / (2 * S * E_w - 0.2 * P)
end

"""
    pressure_vessel_thicknesses(P_bar::Real, D_i_m::Real;
                                 S_MPa::Real=137.0, E_w::Real=1.0)

Compute both cylinder and hemisphere wall thicknesses for a pressure vessel.

# Arguments
- `P_bar`: Internal design pressure in bar
- `D_i_m`: Inner diameter in meters
- `S_MPa`: Maximum allowable stress in MPa (default: 137 MPa for stainless steel 316L)
- `E_w`: Joint efficiency factor (default: 1.0 for full radiographic examination)

# Returns
- Named tuple (t_cylinder_mm, t_hemisphere_mm)
"""
function pressure_vessel_thicknesses(P_bar::Real, D_i_m::Real;
                                      S_MPa::Real=137.0, E_w::Real=1.0)
    # Convert units
    P_MPa = P_bar * 0.1  # bar → MPa
    R_i_mm = D_i_m * 1000.0 / 2.0  # m → mm, diameter → radius

    t_cyl = cylinder_wall_thickness(P_MPa, R_i_mm, S_MPa, E_w)
    t_hem = hemisphere_wall_thickness(P_MPa, R_i_mm, S_MPa, E_w)

    return (t_cylinder_mm=t_cyl, t_hemisphere_mm=t_hem)
end

# =============================================================================
# ITACA Full Dimensions
# =============================================================================

# Component thicknesses (cm)
# Radial components
const BFD_thickness_cm = 1.0        # Buffer Field region
const clearance_cm = 1.0            # Clearances
const ICS_thickness_cm = 15.0       # Inner Copper Shield

# Axial components
const DSP_LG_FATGEM_cm = 10.0       # DSP + Light Guide + FAT-GEM combined
const MRS_thickness_cm = 10.0       # MRS area (below cathode)

# Titanium Grade 5 (Ti-6Al-4V) allowable stress
const S_Ti_Grade5_MPa = 240.0       # MPa (typical allowable stress)

"""
    itaca_full_dimensions(; D_fid_m::Real=D_fid_default,
                           L_TPC_m::Real=L_fid_default,
                           P_bar::Real=15.0,
                           S_MPa::Real=S_Ti_Grade5_MPa,
                           E_w::Real=1.0)

Compute the full dimensions of the ITACA detector including all components.

## Radial structure (from center outward):
1. TPC (fiducial radius R_fid = D_fid/2)
2. BFD (Buffer Field region): 1 cm
3. Clearances: 1 cm
4. ICS (Inner Copper Shield): 15 cm
5. PV (Pressure Vessel): Titanium Grade 5, thickness computed via ASME

## Axial structure (Z, from top to bottom):
1. ICS (top): 15 cm
2. DSP + Light Guide + FAT-GEM: 10 cm (combined)
3. TPC: L_TPC (drift length, 150 cm)
4. Cathode (negligible thickness)
5. MRS area: 10 cm
6. ICS (bottom): 15 cm

# Arguments
- `D_fid_m`: Fiducial diameter in meters (default: 3.2 m)
- `L_TPC_m`: TPC length in meters (default: 1.5 m)
- `P_bar`: Operating pressure in bar (default: 15)
- `S_MPa`: Material allowable stress in MPa (default: 240 MPa for Ti Grade 5)
- `E_w`: Joint efficiency factor (default: 1.0)

# Returns
- Named tuple with all dimensions
"""
function itaca_full_dimensions(; D_fid_m::Real=D_fid_default,
                                 L_TPC_m::Real=L_fid_default,
                                 P_bar::Real=15.0,
                                 S_MPa::Real=S_Ti_Grade5_MPa,
                                 E_w::Real=1.0)
    # =========================================================================
    # Radial dimensions
    # =========================================================================
    R_fid_cm = D_fid_m * 100.0 / 2.0  # Fiducial radius in cm

    # Add components radially: TPC → BFD → Clearances → ICS → PV
    R_BFD_cm = R_fid_cm + BFD_thickness_cm
    R_clearance_cm = R_BFD_cm + clearance_cm
    R_ICS_cm = R_clearance_cm + ICS_thickness_cm

    # Inner radius of pressure vessel = outer radius of ICS
    R_PV_inner_cm = R_ICS_cm
    D_PV_inner_m = 2.0 * R_PV_inner_cm / 100.0  # Convert to meters

    # Compute PV wall thickness using ASME formula
    thicknesses = pressure_vessel_thicknesses(P_bar, D_PV_inner_m; S_MPa=S_MPa, E_w=E_w)
    t_cyl_mm = thicknesses.t_cylinder_mm
    t_hem_mm = thicknesses.t_hemisphere_mm

    # Outer radius of PV
    R_PV_outer_cm = R_PV_inner_cm + t_cyl_mm / 10.0  # mm to cm
    D_PV_outer_m = 2.0 * R_PV_outer_cm / 100.0

    # =========================================================================
    # Axial dimensions (Z direction)
    # =========================================================================
    L_TPC_cm = L_TPC_m * 100.0  # TPC length in cm

    # Total internal length (from ICS top to ICS bottom)
    L_internal_cm = (ICS_thickness_cm +      # Top ICS
                     DSP_LG_FATGEM_cm +       # DSP + LG + FAT-GEM
                     L_TPC_cm +               # TPC drift region
                     MRS_thickness_cm +       # MRS area
                     ICS_thickness_cm)        # Bottom ICS

    L_internal_m = L_internal_cm / 100.0

    # Total external length (add hemisphere caps)
    L_external_m = L_internal_m + 2.0 * (t_hem_mm / 1000.0)

    return (
        # Fiducial
        D_fid_m = D_fid_m,
        R_fid_cm = R_fid_cm,
        L_TPC_m = L_TPC_m,
        L_TPC_cm = L_TPC_cm,

        # Radial build-up (cm)
        R_BFD_cm = R_BFD_cm,
        R_clearance_cm = R_clearance_cm,
        R_ICS_cm = R_ICS_cm,
        R_PV_inner_cm = R_PV_inner_cm,
        R_PV_outer_cm = R_PV_outer_cm,

        # Diameters (m)
        D_PV_inner_m = D_PV_inner_m,
        D_PV_outer_m = D_PV_outer_m,

        # PV wall thicknesses (mm)
        t_cyl_mm = t_cyl_mm,
        t_hem_mm = t_hem_mm,

        # Axial dimensions
        L_internal_cm = L_internal_cm,
        L_internal_m = L_internal_m,
        L_external_m = L_external_m,

        # Component thicknesses (for reference)
        BFD_cm = BFD_thickness_cm,
        clearance_cm = clearance_cm,
        ICS_cm = ICS_thickness_cm,
        DSP_LG_FATGEM_cm = DSP_LG_FATGEM_cm,
        MRS_cm = MRS_thickness_cm,

        # Material and pressure
        P_bar = P_bar,
        S_MPa = S_MPa,
        E_w = E_w
    )
end

"""
    print_itaca_dimensions(dims::NamedTuple)

Print a formatted summary of ITACA detector dimensions.
"""
function print_itaca_dimensions(dims::NamedTuple)
    println("\n" * "=" ^ 60)
    println("          ITACA Full Detector Dimensions")
    println("=" ^ 60)

    println("\n--- Design Parameters ---")
    @printf("  Operating pressure:    P = %.0f bar\n", dims.P_bar)
    @printf("  PV Material:           Ti Grade 5 (S = %.0f MPa)\n", dims.S_MPa)
    @printf("  Joint efficiency:      E_w = %.1f\n", dims.E_w)

    println("\n--- Fiducial Volume ---")
    @printf("  Fiducial diameter:     D_fid = %.2f m\n", dims.D_fid_m)
    @printf("  Fiducial radius:       R_fid = %.1f cm\n", dims.R_fid_cm)
    @printf("  TPC length (drift):    L_TPC = %.2f m (%.0f cm)\n", dims.L_TPC_m, dims.L_TPC_cm)

    println("\n--- Radial Structure (from center outward) ---")
    @printf("  1. TPC (fiducial):     R = %.1f cm\n", dims.R_fid_cm)
    @printf("  2. + BFD (%.0f cm):       R = %.1f cm\n", dims.BFD_cm, dims.R_BFD_cm)
    @printf("  3. + Clearances (%.0f cm): R = %.1f cm\n", dims.clearance_cm, dims.R_clearance_cm)
    @printf("  4. + ICS (%.0f cm):      R = %.1f cm  ← PV inner radius\n", dims.ICS_cm, dims.R_ICS_cm)
    @printf("  5. + PV wall:          R = %.1f cm  (t = %.2f mm)\n", dims.R_PV_outer_cm, dims.t_cyl_mm)

    println("\n--- Axial Structure (Z, top to bottom) ---")
    @printf("  1. ICS (top):          %.0f cm\n", dims.ICS_cm)
    @printf("  2. DSP+LG+FAT-GEM:     %.0f cm\n", dims.DSP_LG_FATGEM_cm)
    @printf("  3. TPC (drift):        %.0f cm\n", dims.L_TPC_cm)
    @printf("  4. Cathode:            ~ 0 cm\n")
    @printf("  5. MRS area:           %.0f cm\n", dims.MRS_cm)
    @printf("  6. ICS (bottom):       %.0f cm\n", dims.ICS_cm)
    @printf("  ---------------------------------\n")
    @printf("     Total internal:     %.0f cm = %.2f m\n", dims.L_internal_cm, dims.L_internal_m)

    println("\n--- Pressure Vessel Dimensions ---")
    @printf("  Inner diameter:        D_i = %.3f m (%.1f cm)\n", dims.D_PV_inner_m, dims.R_PV_inner_cm * 2)
    @printf("  Outer diameter:        D_o = %.3f m (%.1f cm)\n", dims.D_PV_outer_m, dims.R_PV_outer_cm * 2)
    @printf("  Cylinder wall:         t_cyl = %.2f mm\n", dims.t_cyl_mm)
    @printf("  Hemisphere wall:       t_hem = %.2f mm\n", dims.t_hem_mm)
    @printf("  Internal length:       L_i = %.2f m\n", dims.L_internal_m)
    @printf("  External length:       L_o = %.2f m (with caps)\n", dims.L_external_m)
end

# =============================================================================
# Ion Diffusion Functions
# =============================================================================

"""
    transverse_diffusion(L::Real, E::Real; T::Real=300.0)

Compute the transverse diffusion σ of thermal ions drifting a distance L
in an electric field E.

For thermal ions, the Einstein relation gives D = μ k_B T / q_e.
The transverse spread after drifting distance L is:

    σ = √(2 D t) = √(2 k_B T L / (q_e E))

Note: This assumes ions remain thermal (no field heating).

# Arguments
- `L`: Drift length in cm
- `E`: Electric field in V/cm
- `T`: Temperature in K (default: 300 K)

# Returns
- Transverse diffusion σ in cm
"""
function transverse_diffusion(L::Real, E::Real; T::Real=300.0)
    # Convert E from V/cm to V/m for SI units
    E_SI = E * 100.0  # V/m
    # Convert L from cm to m
    L_SI = L / 100.0  # m

    σ_SI = sqrt(2 * k_B * T * L_SI / (q_e * E_SI))
    return σ_SI * 100.0  # Convert back to cm
end

"""
    transverse_diffusion_mm(L_mm::Real, E::Real; T::Real=300.0)

Same as `transverse_diffusion` but with L in mm and returns σ in mm.

# Arguments
- `L_mm`: Drift length in mm
- `E`: Electric field in V/cm
- `T`: Temperature in K (default: 300 K)

# Returns
- Transverse diffusion σ in mm
"""
function transverse_diffusion_mm(L_mm::Real, E::Real; T::Real=300.0)
    L_cm = L_mm / 10.0
    σ_cm = transverse_diffusion(L_cm, E; T=T)
    return σ_cm * 10.0  # cm → mm
end

# =============================================================================
# FAT-GEM Electrode Voltage Functions
# =============================================================================

"""
    voltage_electrode_1(V_c::Real, E::Real, L_m::Real)

Compute the voltage at electrode 1 (entry electrode of FAT-GEM).

    V_electrode_1 = V_c + E × L

where E is in V/cm and L is the drift length.

# Arguments
- `V_c`: Cathode voltage in kV
- `E`: Electric field in V/cm
- `L_m`: Drift length in meters

# Returns
- Voltage at electrode 1 in kV
"""
function voltage_electrode_1(V_c::Real, E::Real, L_m::Real)
    L_cm = L_m * 100.0  # m → cm
    V_drift = E * L_cm / 1000.0  # V → kV
    return V_c + V_drift
end

"""
    voltage_electrode_2(V_c::Real, E::Real, L_m::Real; V_fatgem::Real=V_FATGEM)

Compute the voltage at electrode 2 (exit electrode of FAT-GEM).

    V_electrode_2 = V_electrode_1 + V_FATGEM

# Arguments
- `V_c`: Cathode voltage in kV
- `E`: Electric field in V/cm
- `L_m`: Drift length in meters
- `V_fatgem`: Voltage across FAT-GEM in kV (default: 15 kV)

# Returns
- Voltage at electrode 2 in kV
"""
function voltage_electrode_2(V_c::Real, E::Real, L_m::Real; V_fatgem::Real=V_FATGEM)
    return voltage_electrode_1(V_c, E, L_m) + V_fatgem
end


"""
MARS Kinematics Module
======================
Computes rotation time, settling time, dead zone, and fiducial efficiency
for the MARS positioning system in the ITACA detector.

All SI units unless noted.
"""

# ─── Moment of inertia ────────────────────────────────────────────────────

"""
    I_arms(μ, R)

Moment of inertia of two uniform arms, each of linear mass density μ (kg/m)
and length R (m), pivoting at one end:

    I_arm = (1/3) μ R³
    I_arms = (2/3) μ R³
"""
I_arms(μ, R) = (2/3) * μ * R^3

"""
    I_plates(m_plate, R)

Moment of inertia of two point-mass ion plates, each of mass m_plate (kg),
at worst-case radial distance R (m):

    I_plates = 2 m_plate R²
"""
I_plates(m_plate, R) = 2 * m_plate * R^2

"""
    I_total(μ, m_plate, R)

Total moment of inertia of the rotating assembly:

    I_total = I_arms + I_plates = (2/3) μ R³ + 2 m_plate R²
"""
I_total(μ, m_plate, R) = I_arms(μ, R) + I_plates(m_plate, R)

# ─── Inertia and drag terms ───────────────────────────────────────────────

"""
    inertia_term(μ, m_plate, R; Δθ=π/2)

    𝓘 = 4 I_total Δθ
"""
inertia_term(μ, m_plate, R; Δθ=π/2) = 4 * I_total(μ, m_plate, R) * Δθ

"""
    drag_term(R; ρ=87.0, Cd=0.1, d_wake=9.6e-3, Δθ=π/2)

    𝓓 = ρ Cd d_wake R⁴ Δθ²
"""
drag_term(R; ρ=87.0, Cd=0.1, d_wake=9.6e-3, Δθ=π/2) = ρ * Cd * d_wake * R^4 * Δθ^2

# ─── Rotation time ────────────────────────────────────────────────────────

"""
    rotation_time(μ, m_plate, R; τ_motor=140.0, ρ=87.0, Cd=0.1, d_wake=9.6e-3, Δθ=π/2)

    t_rot = √((𝓘 + 𝓓) / τ_motor)
"""
function rotation_time(μ, m_plate, R;
                       τ_motor=140.0, ρ=87.0, Cd=0.1, d_wake=9.6e-3, Δθ=π/2)
    𝓘 = inertia_term(μ, m_plate, R; Δθ)
    𝓓 = drag_term(R; ρ, Cd, d_wake, Δθ)
    sqrt((𝓘 + 𝓓) / τ_motor)
end

# ─── Wake velocity and settling time ──────────────────────────────────────

"""
    wake_velocity(R; d_wake=9.6e-3, Δθ=π/2, τ_motor=140.0, μ=1.62, m_plate=0.6)

    α = 4 Δθ / t_rot²
    u₀ = √(α · R · d_wake)
"""
function wake_velocity(R; d_wake=9.6e-3, Δθ=π/2, τ_motor=140.0, μ=1.62, m_plate=0.6)
    t_rot = rotation_time(μ, m_plate, R; τ_motor)
    α = 4Δθ / t_rot^2
    sqrt(α * R * d_wake)
end

"""
    settling_time(R; v_d=0.10, n=1.2, d_wake=9.6e-3, Δθ=π/2, τ_motor=140.0, μ=1.62, m_plate=0.6)

    τ₀ = d_wake / u₀
    t_settle = τ₀ × [(u₀/v_d)^(2/n) − 1]
"""
function settling_time(R; v_d=0.10, n=1.2, d_wake=9.6e-3,
                       Δθ=π/2, τ_motor=140.0, μ=1.62, m_plate=0.6)
    u₀ = wake_velocity(R; d_wake, Δθ, τ_motor, μ, m_plate)
    τ₀ = d_wake / u₀
    τ₀ * ((u₀ / v_d)^(2/n) - 1)
end

# ─── Dead zone and fiducial efficiency ────────────────────────────────────

"""
    dead_zone(μ, m_plate, R; v_d=0.10, ...)

    Z_dead = v_d × (t_rot + t_settle)
"""
function dead_zone(μ, m_plate, R; v_d=0.10, n=1.2, d_wake=9.6e-3,
                   Δθ=π/2, τ_motor=140.0, ρ=87.0, Cd=0.1)
    t_rot = rotation_time(μ, m_plate, R; τ_motor, ρ, Cd, d_wake, Δθ)
    t_set = settling_time(R; v_d, n, d_wake, Δθ, τ_motor, μ, m_plate)
    v_d * (t_rot + t_set)
end

"""
    fiducial_efficiency(μ, m_plate, R, L_fid; v_d=0.10, ...)

    ε_geo = (L_fid − Z_dead) / L_fid
"""
function fiducial_efficiency(μ, m_plate, R, L_fid; v_d=0.10, n=1.2, d_wake=9.6e-3,
                             Δθ=π/2, τ_motor=140.0, ρ=87.0, Cd=0.1)
    Z = dead_zone(μ, m_plate, R; v_d, n, d_wake, Δθ, τ_motor, ρ, Cd)
    (L_fid - Z) / L_fid
end

# =============================================================================
# Dense Silicon Plane (DSP) Functions
# =============================================================================

# DSP default parameters
const DSP_pitch_mm = 10.0       # mm (pitch between SiPM centers)
const SiPM_size_mm = 6.0        # mm (SiPM active area is 6×6 mm²)

"""
    dsp_n_sipms_1d(D_cm::Real; pitch_mm::Real=DSP_pitch_mm)

Compute the number of SiPMs along one dimension (diameter) of the DSP.

For a circular DSP of diameter D, SiPMs are arranged on a square grid with pitch p.
The number along one dimension is approximately D/p.

# Arguments
- `D_cm`: DSP diameter in cm
- `pitch_mm`: SiPM pitch in mm (default: 10 mm)

# Returns
- Number of SiPMs along one dimension (integer)
"""
function dsp_n_sipms_1d(D_cm::Real; pitch_mm::Real=DSP_pitch_mm)
    D_mm = D_cm * 10.0
    return floor(Int, D_mm / pitch_mm)
end

"""
    dsp_n_channels(D_cm::Real; pitch_mm::Real=DSP_pitch_mm)

Compute the total number of SiPM channels in the circular DSP.

SiPMs are arranged on a square grid. Only those whose centers fall within
the circular DSP area are counted.

# Arguments
- `D_cm`: DSP diameter in cm
- `pitch_mm`: SiPM pitch in mm (default: 10 mm)

# Returns
- Total number of SiPM channels (integer)
"""
function dsp_n_channels(D_cm::Real; pitch_mm::Real=DSP_pitch_mm)
    D_mm = D_cm * 10.0
    R_mm = D_mm / 2.0

    # Number of grid points along each axis (centered at 0)
    n_half = floor(Int, R_mm / pitch_mm)

    # Count SiPMs whose centers are within the circular area
    count = 0
    for i in -n_half:n_half
        for j in -n_half:n_half
            x = i * pitch_mm
            y = j * pitch_mm
            if x^2 + y^2 <= R_mm^2
                count += 1
            end
        end
    end
    return count
end

"""
    dsp_coverage(; pitch_mm::Real=DSP_pitch_mm, sipm_size_mm::Real=SiPM_size_mm)

Compute the geometric coverage (fill factor) of SiPMs on the DSP.

Coverage = (SiPM area) / (pitch²) = (sipm_size)² / (pitch)²

# Arguments
- `pitch_mm`: SiPM pitch in mm (default: 10 mm)
- `sipm_size_mm`: SiPM active area side length in mm (default: 6 mm)

# Returns
- Coverage as a fraction (0 to 1)
"""
function dsp_coverage(; pitch_mm::Real=DSP_pitch_mm, sipm_size_mm::Real=SiPM_size_mm)
    return (sipm_size_mm / pitch_mm)^2
end

"""
    dsp_total_sipm_area_m2(D_cm::Real; pitch_mm::Real=DSP_pitch_mm, sipm_size_mm::Real=SiPM_size_mm)

Compute the total active SiPM area in the DSP.

# Arguments
- `D_cm`: DSP diameter in cm
- `pitch_mm`: SiPM pitch in mm (default: 10 mm)
- `sipm_size_mm`: SiPM active area side length in mm (default: 6 mm)

# Returns
- Total SiPM active area in m²
"""
function dsp_total_sipm_area_m2(D_cm::Real; pitch_mm::Real=DSP_pitch_mm, sipm_size_mm::Real=SiPM_size_mm)
    n_channels = dsp_n_channels(D_cm; pitch_mm=pitch_mm)
    sipm_area_mm2 = sipm_size_mm^2
    return n_channels * sipm_area_mm2 * 1e-6  # mm² to m²
end

"""
    print_dsp_summary(D_cm::Real; pitch_mm::Real=DSP_pitch_mm, sipm_size_mm::Real=SiPM_size_mm)

Print a summary of DSP parameters and computed values.

# Arguments
- `D_cm`: DSP diameter in cm
- `pitch_mm`: SiPM pitch in mm (default: 10 mm)
- `sipm_size_mm`: SiPM active area side length in mm (default: 6 mm)
"""
function print_dsp_summary(D_cm::Real; pitch_mm::Real=DSP_pitch_mm, sipm_size_mm::Real=SiPM_size_mm)
    n_channels = dsp_n_channels(D_cm; pitch_mm=pitch_mm)
    coverage = dsp_coverage(; pitch_mm=pitch_mm, sipm_size_mm=sipm_size_mm)
    total_area = dsp_total_sipm_area_m2(D_cm; pitch_mm=pitch_mm, sipm_size_mm=sipm_size_mm)
    dsp_area = π * (D_cm / 200.0)^2  # m²

    println("\n" * "=" ^ 50)
    println("     Dense Silicon Plane (DSP) Summary")
    println("=" ^ 50)
    @printf("\n--- Input Parameters ---\n")
    @printf("  DSP diameter:       D = %.0f cm (%.2f m)\n", D_cm, D_cm/100)
    @printf("  SiPM pitch:         p = %.1f mm\n", pitch_mm)
    @printf("  SiPM size:          %.1f × %.1f mm²\n", sipm_size_mm, sipm_size_mm)

    @printf("\n--- Computed Values ---\n")
    @printf("  Number of channels: N = %d\n", n_channels)
    @printf("  Coverage:           %.1f %%\n", coverage * 100)
    @printf("  DSP total area:     %.2f m²\n", dsp_area)
    @printf("  SiPM active area:   %.2f m²\n", total_area)
end

# =============================================================================
# Plotting Functions
# =============================================================================

# =============================================================================
# Electron Diffusion in Pure Xenon
# =============================================================================

# Transverse diffusion coefficient for electrons in pure Xe
# σ_T = σ_coeff × √(L / P), with σ_coeff in mm × √(bar/cm)
const σ_T_electron_coefficient = 3.5  # mm × √(bar/cm)

"""
    electron_transverse_diffusion_mm(L_cm::Real, P_bar::Real; σ_coeff::Real=σ_T_electron_coefficient)

Compute transverse diffusion of electrons in pure xenon.

For electrons in pure Xe, diffusion scales as:
    σ_T = σ_coeff × √(L / P)

where σ_coeff = 3.5 mm × √(bar/cm).

Dimensional analysis:
    [mm × √(bar/cm)] × √([cm] / [bar]) = mm × √(bar × cm / (cm × bar)) = mm ✓

# Arguments
- `L_cm`: Drift length in cm
- `P_bar`: Pressure in bar
- `σ_coeff`: Diffusion coefficient in mm × √(bar/cm) (default: 3.5)

# Returns
- Transverse diffusion σ in mm
"""
function electron_transverse_diffusion_mm(L_cm::Real, P_bar::Real; σ_coeff::Real=σ_T_electron_coefficient)
    return σ_coeff * sqrt(L_cm / P_bar)
end

"""
    plot_electron_ion_diffusion(; E::Real=300.0, P_bar::Real=15.0,
                                  T::Real=300.0, L_max_cm::Real=150.0)

Create a plot comparing electron and ion transverse diffusion vs drift length.

This plot demonstrates the anti-correlation: electrons have much larger diffusion
(~10× at full drift) compared to thermal ions.

# Arguments
- `E`: Electric field in V/cm (default: 300)
- `P_bar`: Pressure in bar (default: 15, used for labeling)
- `T`: Temperature in K (default: 300)
- `L_max_cm`: Maximum drift length in cm (default: 150)

# Returns
- Plot object
"""
function plot_electron_ion_diffusion(; E::Real=300.0, P_bar::Real=15.0,
                                       T::Real=300.0, L_max_cm::Real=150.0)
    L_range_cm = range(1.0, L_max_cm, length=100)  # Start from 1 cm to avoid sqrt(0)

    # Ion diffusion (thermal limit)
    σ_ion = [transverse_diffusion_mm(L * 10.0, E; T=T) for L in L_range_cm]  # L in mm

    # Electron diffusion (pure Xe) - pass pressure
    σ_electron = [electron_transverse_diffusion_mm(L, P_bar) for L in L_range_cm]

    p = plot(
        xlabel = "Drift length L (cm)",
        ylabel = "Transverse diffusion σ_T (mm)",
        title = "Electron vs Ion Diffusion (E = $(Int(E)) V/cm, P = $(Int(P_bar)) bar)",
        legend = :topleft,
        grid = true,
        size = (800, 550),
        ylims = (0, maximum(σ_electron) * 1.1)
    )

    # Plot electron diffusion (dashed, red)
    plot!(p, L_range_cm, σ_electron,
        label = "Electrons (pure Xe)",
        linewidth = 2.5,
        linestyle = :dash,
        linecolor = :red
    )

    # Plot ion diffusion (solid, blue)
    plot!(p, L_range_cm, σ_ion,
        label = "Ions (thermal)",
        linewidth = 2.5,
        linestyle = :solid,
        linecolor = :blue
    )

    return p
end

"""
    plot_ion_diffusion_vs_length(; E_values::Vector{<:Real}=[200.0, 300.0, 400.0],
                                   T::Real=300.0,
                                   L_min::Real=0.0, L_max::Real=1.5)

Create a plot of ion transverse diffusion σ (mm) vs drift length L (m)
for multiple electric field values.

# Arguments
- `E_values`: Vector of electric field values in V/cm (default: [200, 300, 400])
- `T`: Temperature in K (default: 300)
- `L_min`: Minimum drift length in meters (default: 0)
- `L_max`: Maximum drift length in meters (default: 1.5)

# Returns
- Plot object
"""
function plot_ion_diffusion_vs_length(; E_values::Vector{<:Real}=[200.0, 300.0, 400.0],
                                        T::Real=300.0,
                                        L_min::Real=0.0, L_max::Real=1.5)
    L_range = range(L_min, L_max, length=100)
    colors = [:blue, :red, :green]

    p = plot(
        xlabel = "Drift length L (m)",
        ylabel = "Transverse diffusion σ (mm)",
        title = "Ion Transverse Diffusion (T = $T K)",
        legend = :topleft,
        grid = true,
        size = (800, 500)
    )

    for (i, E) in enumerate(E_values)
        # Convert L from m to mm for the diffusion function
        σ_values = [transverse_diffusion_mm(L * 1000.0, E; T=T) for L in L_range]
        plot!(p, L_range, σ_values,
            label = "E = $(Int(E)) V/cm",
            linewidth = 2.5,
            linecolor = colors[mod1(i, length(colors))]
        )
    end

    return p
end

"""
    plot_electrode_voltage_vs_length(; E_values::Vector{<:Real}=[200.0, 300.0, 400.0],
                                       V_c::Real=0.0,
                                       L_min::Real=0.0, L_max::Real=1.5)

Create a plot of electrode 1 voltage (kV) vs drift length L (m)
for multiple electric field values.

# Arguments
- `E_values`: Vector of electric field values in V/cm (default: [200, 300, 400])
- `V_c`: Cathode voltage in kV (default: 0)
- `L_min`: Minimum drift length in meters (default: 0)
- `L_max`: Maximum drift length in meters (default: 1.5)

# Returns
- Plot object
"""
function plot_electrode_voltage_vs_length(; E_values::Vector{<:Real}=[200.0, 300.0, 400.0],
                                            V_c::Real=0.0,
                                            L_min::Real=0.0, L_max::Real=1.5)
    L_range = range(L_min, L_max, length=100)
    colors = [:blue, :red, :green]

    p = plot(
        xlabel = "Drift length L (m)",
        ylabel = "Electrode 1 voltage (kV)",
        title = "FAT-GEM Entry Electrode Voltage (V_c = $V_c kV)",
        legend = :topleft,
        grid = true,
        size = (800, 500)
    )

    for (i, E) in enumerate(E_values)
        V_values = [voltage_electrode_1(V_c, E, L) for L in L_range]
        plot!(p, L_range, V_values,
            label = "E = $(Int(E)) V/cm",
            linewidth = 2.5,
            linecolor = colors[mod1(i, length(colors))]
        )
    end

    return p
end

# =============================================================================
# Fiducial Mass Functions
# =============================================================================

"""
    fiducial_mass_kg(D_fid::Real, L_fid::Real; ρ::Real=ρ_Xe_15bar)

Compute the fiducial mass of the ITACA detector.

The fiducial volume is a cylinder with diameter D_fid and length L_fid.
    V = π × (D/2)² × L

# Arguments
- `D_fid`: Fiducial diameter in meters
- `L_fid`: Fiducial length in meters
- `ρ`: Xenon density in kg/m³ (default: 87 kg/m³ at 15 bar)

# Returns
- Fiducial mass in kg
"""
function fiducial_mass_kg(D_fid::Real, L_fid::Real; ρ::Real=ρ_Xe_15bar)
    r_fid = D_fid / 2.0
    V_fid = π * r_fid^2 * L_fid  # m³
    return ρ * V_fid
end

"""
    fiducial_mass_ton(D_fid::Real, L_fid::Real; ρ::Real=ρ_Xe_15bar)

Same as `fiducial_mass_kg` but returns mass in metric tons.
"""
function fiducial_mass_ton(D_fid::Real, L_fid::Real; ρ::Real=ρ_Xe_15bar)
    return fiducial_mass_kg(D_fid, L_fid; ρ=ρ) / 1000.0
end

"""
    plot_fiducial_mass_vs_diameter(; L_values::Vector{<:Real}=[1.0, 1.5, 2.0],
                                     D_min::Real=1.0, D_max::Real=4.0,
                                     ρ::Real=ρ_Xe_15bar)

Create a plot of fiducial mass (ton) vs fiducial diameter (m)
for multiple fiducial lengths.

# Arguments
- `L_values`: Vector of fiducial lengths in meters (default: [1.0, 1.5, 2.0])
- `D_min`: Minimum fiducial diameter in meters (default: 1.0)
- `D_max`: Maximum fiducial diameter in meters (default: 4.0)
- `ρ`: Xenon density in kg/m³ (default: 87 kg/m³ at 15 bar)

# Returns
- Plot object
"""
function plot_fiducial_mass_vs_diameter(; L_values::Vector{<:Real}=[1.0, 1.5, 2.0],
                                          D_min::Real=1.0, D_max::Real=4.0,
                                          ρ::Real=ρ_Xe_15bar)
    D_range = range(D_min, D_max, length=100)
    colors = [:blue, :red, :green, :orange, :purple]

    p = plot(
        xlabel = "Fiducial diameter (m)",
        ylabel = "Fiducial mass (ton)",
        title = "ITACA Fiducial Mass (ρ_Xe = $ρ kg/m³ at 15 bar)",
        legend = :topleft,
        grid = true,
        size = (800, 500)
    )

    for (i, L) in enumerate(L_values)
        mass_values = [fiducial_mass_ton(D, L; ρ=ρ) for D in D_range]
        plot!(p, D_range, mass_values,
            label = "L_fid = $L m",
            linewidth = 2.5,
            linecolor = colors[mod1(i, length(colors))]
        )
    end

    return p
end


"""
Arm structural parameters
==========================
Parameterized in terms of number of longerons N and tube cross-section A_tube.
"""

# ─── Material properties ──────────────────────────────────────────────────

const ρ_Ti = 4430.0   # kg/m³ — Ti Gr.5 density

# ─── Tube geometry ────────────────────────────────────────────────────────

"""
    A_tube(d_o, t_w)

Cross-sectional area of a single tube (m²).
    d_o : outer diameter (m)
    t_w : wall thickness (m)
"""
A_tube(d_o, t_w) = (π/4) * (d_o^2 - (d_o - 2t_w)^2)

"""
    I_tube(d_o, t_w)

Second moment of area of a single tube about its own centroid (m⁴).
"""
I_tube(d_o, t_w) = (π/64) * (d_o^4 - (d_o - 2t_w)^4)

# ─── Linear mass density ─────────────────────────────────────────────────

"""
    μ_arm(N, d_o, t_w; ρ=ρ_Ti)

Linear mass density of the arm (kg/m), dominated by N Ti longerons.

    μ = ρ × N × A_tube
"""
μ_arm(N, d_o, t_w; ρ=ρ_Ti) = ρ * N * A_tube(d_o, t_w)

# ─── Arm mass ─────────────────────────────────────────────────────────────

"""
    m_arm(N, d_o, t_w, R; ρ=ρ_Ti)

Total mass of one arm (kg).

    m = μ × R
"""
m_arm(N, d_o, t_w, R; ρ=ρ_Ti) = μ_arm(N, d_o, t_w; ρ) * R

# ─── Section modulus ──────────────────────────────────────────────────────

"""
    I_section(N, d_o, t_w, positions)

Second moment of area of the full longeron truss (m⁴).
Uses parallel axis theorem.

    positions : vector of N vertical offsets from neutral axis (m)
    I_section = N × I_tube + A_tube × Σ yₖ²
"""
function I_section(N, d_o, t_w, positions)
    It = I_tube(d_o, t_w)
    At = A_tube(d_o, t_w)
    N * It + At * sum(y^2 for y in positions)
end

"""
    Z_section(N, d_o, t_w, positions)

Section modulus (m³) to the outermost longeron.

    Z = I_section / y_max
"""
function Z_section(N, d_o, t_w, positions)
    I = I_section(N, d_o, t_w, positions)
    y_max = maximum(abs.(positions))
    I / y_max
end


"""
Structural analysis  - Load cases
=================================
Simply supported beam: pin at hub (r=0), roller at Cu rail (r=R).
"""

const g = 9.81  # m/s²

# ─── Load Case 1: Tangential bending (during rotation) ────────────────────

"""
    M_tang(μ, α, R)

Maximum tangential bending moment (N·m).
Linearly increasing inertial load q(r) = μ α r.
Maximum at r = R/√3.

    M_tang = μ α R³ / (9√3)
"""
M_tang(μ, α, R) = μ * α * R^3 / (9√3)

"""
    r_tang_max(R)

Radial position of maximum tangential moment (m).

    r = R / √3
"""
r_tang_max(R) = R / √3

"""
    M_tang_at(μ, α, R, r)

Tangential bending moment at arbitrary position r (N·m).

    M(r) = (μ α / 6)(R² r − r³)
"""
M_tang_at(μ, α, R, r) = (μ * α / 6) * (R^2 * r - r^3)

# ─── Load Case 2: Gravity bending ────────────────────────────────────────

"""
    M_grav(μ, R)

Maximum gravity bending moment (N·m).
Uniform load q = μ g, maximum at midspan r = R/2.

    M_grav = μ g R² / 8
"""
M_grav(μ, R) = μ * g * R^2 / 8

"""
    M_grav_at(μ, R, r)

Gravity bending moment at arbitrary position r (N·m).

    M(r) = (μ g / 2) r (R − r)
"""
M_grav_at(μ, R, r) = (μ * g / 2) * r * (R - r)

# ─── Stress from bending ─────────────────────────────────────────────────

"""
    σ_bending(M, Z)

Bending stress (Pa).

    σ = M / Z
"""
σ_bending(M, Z) = M / Z


"""
Deflection analysis
====================
Simply supported beam: pin at hub (r=0), roller at Cu rail (r=R).
"""

const E_Ti = 114e9  # Pa — Young's modulus Ti Gr.5

# ─── Flexural rigidity ────────────────────────────────────────────────────

"""
    EI(N, d_o, t_w, positions; E=E_Ti)

Flexural rigidity of the longeron truss (N·m²).

    EI = E × I_section
"""
EI_flex(N, d_o, t_w, positions; E=E_Ti) = E * I_section(N, d_o, t_w, positions)

# ─── Load Case 1: Tangential deflection ──────────────────────────────────

"""
    δ_tang_at(μ, α, R, r, EI)

Tangential deflection at arbitrary position r (m).

    y(r) = (μα / 360EI)(10R²r³ − 3r⁵ − 7R⁴r)
"""
δ_tang_at(μ, α, R, r, EI) = (μ * α / (360 * EI)) * (10R^2 * r^3 - 3r^5 - 7R^4 * r)

"""
    δ_tang_max(μ, α, R, EI)

Maximum tangential deflection (m), at r ≈ 0.519R.

    δ_tang = μ α R⁵ / (153 EI)
"""
δ_tang_max(μ, α, R, EI) = μ * α * R^5 / (153 * EI)

"""
    r_tang_max_deflection(R)

Position of maximum tangential deflection (m).

    r ≈ 0.519 R
"""
r_tang_max_deflection(R) = 0.519 * R

# ─── Load Case 2: Gravity deflection ─────────────────────────────────────

"""
    δ_grav_at(μ, R, r, EI)

Gravity deflection at arbitrary position r (m).

    y(r) = (μg / 24EI) r(R³ − 2Rr² + r³)
"""
δ_grav_at(μ, R, r, EI) = (μ * g / (24 * EI)) * r * (R^3 - 2R * r^2 + r^3)

"""
    δ_grav_max(μ, R, EI)

Maximum gravity deflection at midspan r = R/2 (m).

    δ_grav = 5μgR⁴ / (384 EI)
"""
δ_grav_max(μ, R, EI) = 5 * μ * g * R^4 / (384 * EI)
