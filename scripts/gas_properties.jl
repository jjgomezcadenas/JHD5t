#!/usr/bin/env julia

# Gas properties for Xe/He mixtures
# Computes drift velocity of Xe2+ ions using Blanc's Law

using Printf
using Plots
using StatsPlots
using Statistics

# Mobility values at STP (cm²/V/s) for ions in pure Xe
const μ_Xe2_STP = 0.74    # Xe2+ in pure Xe at STP
const μ_He_STP = 19.4     # Xe2+ in pure He at STP (for mixture calculation)
const μ_Ba2_STP = 0.56    # Ba2+ in pure Xe at STP
const μ_NH4_STP = 1.63    # NH4+ in pure Xe at STP (interpolated from alkaline ion fit)

# Alkaline ion data: ionic radius (Å) and reduced mobility K₀ (cm²/V/s) in Xe
const ALKALINE_IONS = (
    Na = (name = "Na⁺", radius = 1.02, K0 = 2.33),
    K  = (name = "K⁺",  radius = 1.38, K0 = 1.80),
    Rb = (name = "Rb⁺", radius = 1.52, K0 = 1.55),
    Cs = (name = "Cs⁺", radius = 1.67, K0 = 1.34)
)

# NH4+ ionic radius (effectively matches Rb+)
const r_NH4 = 1.48  # Å

# Physical constants
const k_B = 1.380649e-23   # Boltzmann constant (J/K)
const q_e = 1.602176634e-19 # Elementary charge (C)
const R_gas = 8.31446      # Universal gas constant (J/(mol·K))

# Molar masses (g/mol)
const M_Xe = 131.293       # Xenon molar mass
const M_He = 4.0026        # Helium molar mass

# Xenon critical constants (for real gas calculations)
const Tc_Xe = 289.73       # Critical temperature (K)
const Pc_Xe = 58.40e5      # Critical pressure (Pa) = 58.40 bar

# Van der Waals constants for Xenon (calculated from critical constants)
# a = 27 R² Tc² / (64 Pc), b = R Tc / (8 Pc)
const a_Xe = 27 * R_gas^2 * Tc_Xe^2 / (64 * Pc_Xe)  # Pa·m⁶/mol²
const b_Xe = R_gas * Tc_Xe / (8 * Pc_Xe)             # m³/mol

# Helium critical constants
const Tc_He = 5.19         # Critical temperature (K)
const Pc_He = 2.27e5       # Critical pressure (Pa) = 2.27 bar

# Van der Waals constants for Helium
const a_He = 27 * R_gas^2 * Tc_He^2 / (64 * Pc_He)
const b_He = R_gas * Tc_He / (8 * Pc_He)

# Standard conditions
const T_STP = 273.15    # K
const P_STP = 1.0       # bar

"""
    linear_fit(x, y)

Perform a simple linear regression y = a + b*x.

# Returns
- (a, b): intercept and slope
"""
function linear_fit(x, y)
    n = length(x)
    x_mean = mean(x)
    y_mean = mean(y)

    b = sum((x .- x_mean) .* (y .- y_mean)) / sum((x .- x_mean).^2)
    a = y_mean - b * x_mean

    return (a, b)
end

"""
    interpolate_NH4_mobility()

Fit a linear model to alkaline ion mobility vs ionic radius,
then interpolate to get NH4+ mobility.

# Returns
- (μ_NH4, a, b): interpolated mobility, intercept, slope
"""
function interpolate_NH4_mobility()
    radii = [ALKALINE_IONS.Na.radius, ALKALINE_IONS.K.radius,
             ALKALINE_IONS.Rb.radius, ALKALINE_IONS.Cs.radius]
    K0s = [ALKALINE_IONS.Na.K0, ALKALINE_IONS.K.K0,
           ALKALINE_IONS.Rb.K0, ALKALINE_IONS.Cs.K0]

    a, b = linear_fit(radii, K0s)
    μ_NH4 = a + b * r_NH4

    return (μ_NH4, a, b)
end

"""
    plot_alkaline_mobility_fit()

Create a plot of reduced mobility K₀ vs ionic radius for alkaline ions,
with linear fit and NH4+ interpolation.

The plot has two x-axes: one with ionic radius, one with ion names.

# Returns
- Plot object
"""
function plot_alkaline_mobility_fit()
    # Data points
    radii = [ALKALINE_IONS.Na.radius, ALKALINE_IONS.K.radius,
             ALKALINE_IONS.Rb.radius, ALKALINE_IONS.Cs.radius]
    K0s = [ALKALINE_IONS.Na.K0, ALKALINE_IONS.K.K0,
           ALKALINE_IONS.Rb.K0, ALKALINE_IONS.Cs.K0]
    names = [ALKALINE_IONS.Na.name, ALKALINE_IONS.K.name,
             ALKALINE_IONS.Rb.name, ALKALINE_IONS.Cs.name]

    # Fit line
    a, b = linear_fit(radii, K0s)
    r_fit = range(0.9, 1.8, length=100)
    K0_fit = a .+ b .* r_fit

    # NH4+ interpolation point
    K0_NH4 = a + b * r_NH4

    # Create plot
    p = plot(r_fit, K0_fit,
        xlabel = "Ionic radius (Å)",
        ylabel = "Reduced mobility K₀ (cm²/V/s)",
        title = "Alkaline Ion Mobility in Xe",
        label = @sprintf("Linear fit: K₀ = %.2f - %.2f × r", a, abs(b)),
        linewidth = 1.5,
        linecolor = :black,
        linestyle = :dash,
        grid = true,
        legend = :topright,
        size = (800, 500),
        top_margin = 15Plots.mm
    )

    # Plot data points
    scatter!(p, radii, K0s,
        markersize = 8,
        markercolor = :blue,
        label = "Alkaline ions"
    )

    # Plot NH4+ interpolated point
    scatter!(p, [r_NH4], [K0_NH4],
        markersize = 10,
        markercolor = :red,
        markershape = :star5,
        label = @sprintf("NH₄⁺ (interpolated): K₀ = %.2f", K0_NH4)
    )

    # Add ion labels above data points
    for (r, K0, name) in zip(radii, K0s, names)
        annotate!(r, K0 + 0.08, text(name, 9, :center))
    end
    annotate!(r_NH4, K0_NH4 + 0.08, text("NH₄⁺", 9, :center, :red))

    # Add secondary x-axis with ion names at top
    # We'll add text annotations at the top margin
    plot!(p,
        top_margin = 20Plots.mm,
        xlims = (0.9, 1.8)
    )

    # Add ion name labels at top of plot
    y_top = maximum(K0s) + 0.25
    for (r, name) in zip(radii, names)
        annotate!(r, y_top, text(name, 8, :center, :gray))
    end
    annotate!(r_NH4, y_top, text("NH₄⁺", 8, :center, :red))

    return p
end

"""
    mobility_mixture_stp(f_Xe::Real, f_He::Real)

Compute the mobility of Xe2+ ions in a Xe/He mixture at STP using Blanc's Law:

    1/μ_mix = f_Xe/μ_Xe + f_He/μ_He

# Arguments
- `f_Xe`: Fraction of Xenon in the mixture (0 to 1)
- `f_He`: Fraction of Helium in the mixture (0 to 1)

# Returns
- Mobility at STP in cm²/V/s
"""
function mobility_mixture_stp(f_Xe::Real, f_He::Real)
    @assert f_Xe + f_He ≈ 1.0 "Fractions must sum to 1"
    @assert 0 ≤ f_Xe ≤ 1 && 0 ≤ f_He ≤ 1 "Fractions must be between 0 and 1"

    inv_μ = f_Xe / μ_Xe2_STP + f_He / μ_He_STP
    return 1.0 / inv_μ
end

"""
    mobility_at_pressure(μ_stp::Real, P::Real)

Scale mobility from STP to a given pressure (assuming constant temperature).

Mobility scales inversely with pressure: μ(P) = μ_STP × (P_STP / P)

# Arguments
- `μ_stp`: Mobility at STP in cm²/V/s
- `P`: Pressure in bar

# Returns
- Mobility at pressure P in cm²/V/s
"""
function mobility_at_pressure(μ_stp::Real, P::Real)
    return μ_stp * (P_STP / P)
end

"""
    drift_velocity(f_Xe::Real, f_He::Real, P::Real, E::Real)

Compute the drift velocity of Xe2+ ions in a Xe/He mixture.

Uses Blanc's Law to determine mixture mobility at STP, then scales
for pressure and applies the electric field:

    v_d = μ(P) × E

# Arguments
- `f_Xe`: Fraction of Xenon in the mixture (0 to 1)
- `f_He`: Fraction of Helium in the mixture (0 to 1)
- `P`: Pressure in bar
- `E`: Electric field in V/cm

# Returns
- Drift velocity in cm/s
"""
function drift_velocity(f_Xe::Real, f_He::Real, P::Real, E::Real)
    μ_stp = mobility_mixture_stp(f_Xe, f_He)
    μ_P = mobility_at_pressure(μ_stp, P)
    return μ_P * E
end

"""
    drift_velocity_mm_μs(f_Xe::Real, f_He::Real, P::Real, E::Real)

Same as `drift_velocity` but returns result in mm/μs.
"""
function drift_velocity_mm_μs(f_Xe::Real, f_He::Real, P::Real, E::Real)
    v_cm_s = drift_velocity(f_Xe, f_He, P, E)
    return v_cm_s * 1e-5  # cm/s → mm/μs
end



"""
    drift_velocity_ion(μ_stp::Real, P::Real, E::Real)

Compute drift velocity for an ion with given mobility at STP in pure Xe.

# Arguments
- `μ_stp`: Ion mobility at STP in cm²/V/s
- `P`: Pressure in bar
- `E`: Electric field in V/cm

# Returns
- Drift velocity in cm/s
"""
function drift_velocity_ion(μ_stp::Real, P::Real, E::Real)
    μ_P = mobility_at_pressure(μ_stp, P)
    return μ_P * E
end

"""
    drift_velocity_ion_mm_s(μ_stp::Real, P::Real, E::Real)

Same as `drift_velocity_ion` but returns result in mm/s.
"""
function drift_velocity_ion_mm_s(μ_stp::Real, P::Real, E::Real)
    v_cm_s = drift_velocity_ion(μ_stp, P, E)
    return v_cm_s * 10.0  # cm/s → mm/s
end

"""
    plot_drift_velocities_bar(P::Real, E::Real; f_Xe::Real=0.90, f_He::Real=0.10)

Create a bar plot comparing drift velocities of different ions.

Shows drift velocities (in mm/s) for:
- Xe2+ (grouped: pure Xe and Xe/He mixture)
- NH4+ in pure Xe
- Ba2+ in pure Xe

# Arguments
- `P`: Pressure in bar
- `E`: Electric field in V/cm
- `f_Xe`: Xenon fraction for mixture (default: 0.90)
- `f_He`: Helium fraction for mixture (default: 0.10)

# Returns
- Plot object
"""
function plot_drift_velocities_bar(P::Real, E::Real; f_Xe::Real=0.90, f_He::Real=0.10)
    # Calculate drift velocities in mm/s
    v_Xe2_pure = drift_velocity_ion_mm_s(μ_Xe2_STP, P, E)
    v_NH4 = drift_velocity_ion_mm_s(μ_NH4_STP, P, E)
    v_Ba2 = drift_velocity_ion_mm_s(μ_Ba2_STP, P, E)
    v_Xe2_mix = drift_velocity(f_Xe, f_He, P, E) * 10.0  # cm/s → mm/s

    # Use explicit bar positions for precise control
    # Xe2+ gets two bars (positions 1 and 2), NH4+ at 3.5, Ba2+ at 5
    positions = [1, 2, 3.5, 5]
    velocities = [v_Xe2_pure, v_Xe2_mix, v_NH4, v_Ba2]
    colors = [:blue, :red, :blue, :blue]

    p = bar(positions, velocities,
        xlabel = "Ion species",
        ylabel = "Drift velocity (mm/s)",
        title = "Ion Drift Velocities (P = $P bar, E = $E V/cm)",
        bar_width = 0.8,
        fillcolor = colors,
        grid = true,
        legend = :topright,
        size = (600, 400),
        top_margin = 10Plots.mm,
        xticks = ([1.5, 3.5, 5], ["Xe₂⁺", "NH₄⁺", "Ba²⁺"]),
        label = nothing
    )

    # Add legend manually
    bar!([NaN], [NaN], fillcolor = :blue, label = "Pure Xe")
    bar!([NaN], [NaN], fillcolor = :red, label = "Xe/He (90/10)")

    # Add value labels exactly on top of each bar with units
    max_v = maximum(velocities)
    for (pos, v) in zip(positions, velocities)
        annotate!(pos, v + max_v*0.05, text(@sprintf("%.0f mm/s", v), 7, :center))
    end

    return p
end

"""
    drift_time_s(L_mm::Real, μ_stp::Real, P::Real, E::Real)

Compute the drift time for an ion to travel distance L.

# Arguments
- `L_mm`: Drift distance in mm
- `μ_stp`: Ion mobility at STP in cm²/V/s
- `P`: Pressure in bar
- `E`: Electric field in V/cm

# Returns
- Drift time in seconds
"""
function drift_time_s(L_mm::Real, μ_stp::Real, P::Real, E::Real)
    v_mm_s = drift_velocity_ion_mm_s(μ_stp, P, E)
    return L_mm / v_mm_s
end

"""
    plot_drift_time_vs_distance(P::Real, E::Real; L_min::Real=10.0, L_max::Real=1300.0)

Create a plot of drift time (s) vs distance (mm) for each ion species in pure Xe.

# Arguments
- `P`: Pressure in bar
- `E`: Electric field in V/cm
- `L_min`: Minimum drift distance in mm (default: 10)
- `L_max`: Maximum drift distance in mm (default: 1300)

# Returns
- Plot object
"""
function plot_drift_time_vs_distance(P::Real, E::Real; L_min::Real=10.0, L_max::Real=1300.0)
    L_range = range(L_min, L_max, length=100)

    t_Xe2 = [drift_time_s(L, μ_Xe2_STP, P, E) for L in L_range]
    t_NH4 = [drift_time_s(L, μ_NH4_STP, P, E) for L in L_range]
    t_Ba2 = [drift_time_s(L, μ_Ba2_STP, P, E) for L in L_range]

    p = plot(L_range, t_Xe2,
        xlabel = "Drift distance (mm)",
        ylabel = "Drift time (s)",
        title = "Ion Drift Time (P = $P bar, E = $E V/cm)",
        label = "Xe₂⁺",
        linewidth = 2,
        linecolor = :blue,
        grid = true,
        legend = :topleft
    )
    plot!(p, L_range, t_NH4, label = "NH₄⁺", linewidth = 2, linecolor = :green)
    plot!(p, L_range, t_Ba2, label = "Ba²⁺", linewidth = 2, linecolor = :orange)

    return p
end

"""
    plot_arrival_time_difference(P::Real, E::Real; L_min::Real=10.0, L_max::Real=1300.0)

Create a plot of arrival time differences vs drift distance.

Shows how much earlier the faster ions arrive compared to Ba²⁺ (slowest):
- Δt(Ba²⁺ - Xe₂⁺): Xe₂⁺ arrives this much earlier than Ba²⁺
- Δt(Ba²⁺ - NH₄⁺): NH₄⁺ arrives this much earlier than Ba²⁺

Mobility order: NH₄⁺ (1.63) > Xe₂⁺ (0.79) > Ba²⁺ (0.56) cm²/V/s
So NH₄⁺ is fastest, then Xe₂⁺, then Ba²⁺ (slowest).

# Arguments
- `P`: Pressure in bar
- `E`: Electric field in V/cm
- `L_min`: Minimum drift distance in mm (default: 10)
- `L_max`: Maximum drift distance in mm (default: 1300)

# Returns
- Plot object
"""
function plot_arrival_time_difference(P::Real, E::Real; L_min::Real=10.0, L_max::Real=1300.0)
    L_range = range(L_min, L_max, length=100)

    # Time differences: how much earlier faster ions arrive compared to Ba2+ (slowest)
    # Δt = t_Ba2 - t_fast_ion (positive values)
    Δt_Ba2_Xe2 = [drift_time_s(L, μ_Ba2_STP, P, E) - drift_time_s(L, μ_Xe2_STP, P, E) for L in L_range]
    Δt_Ba2_NH4 = [drift_time_s(L, μ_Ba2_STP, P, E) - drift_time_s(L, μ_NH4_STP, P, E) for L in L_range]

    p = plot(L_range, Δt_Ba2_NH4,
        xlabel = "Drift distance (mm)",
        ylabel = "Arrival time difference (s)",
        title = "Ion Δt vs Ba²⁺ (P = $P bar, E = $E V/cm)",
        label = "Δt(Ba²⁺ - NH₄⁺)",
        linewidth = 1.5,
        linecolor = :green,
        grid = true,
        legend = :topleft,
        size = (700, 500),
        ylims = (0, 4)
    )
    plot!(p, L_range, Δt_Ba2_Xe2,
        label = "Δt(Ba²⁺ - Xe₂⁺)",
        linewidth = 1.5,
        linecolor = :blue)

    return p
end

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

# Electron transverse diffusion coefficients D_T (√bar · mm / √cm)
const D_T_Xe_pure = 3.5    # Pure Xe
const D_T_Xe_He = 1.6      # Xe/He (90/10)

"""
    electron_transverse_diffusion_mm(L_mm::Real, P_bar::Real; D_T::Real=D_T_Xe_pure)

Compute the transverse diffusion σ of electrons drifting a distance L
at pressure P.

The formula is: σ_T = D_T × √(L/P)

where D_T is the transverse diffusion coefficient with units √bar · mm / √cm.

# Arguments
- `L_mm`: Drift length in mm
- `P_bar`: Pressure in bar
- `D_T`: Transverse diffusion coefficient in √bar · mm / √cm (default: pure Xe)

# Returns
- Transverse diffusion σ in mm
"""
function electron_transverse_diffusion_mm(L_mm::Real, P_bar::Real; D_T::Real=D_T_Xe_pure)
    L_cm = L_mm / 10.0  # Convert mm to cm
    σ_mm = D_T * sqrt(L_cm / P_bar)
    return σ_mm
end

"""
    electron_transverse_diffusion_pure_Xe_mm(L_mm::Real, P_bar::Real)

Compute electron transverse diffusion in pure Xe.

# Arguments
- `L_mm`: Drift length in mm
- `P_bar`: Pressure in bar

# Returns
- Transverse diffusion σ in mm
"""
function electron_transverse_diffusion_pure_Xe_mm(L_mm::Real, P_bar::Real)
    return electron_transverse_diffusion_mm(L_mm, P_bar; D_T=D_T_Xe_pure)
end

"""
    electron_transverse_diffusion_Xe_He_mm(L_mm::Real, P_bar::Real)

Compute electron transverse diffusion in Xe/He (90/10) mixture.

# Arguments
- `L_mm`: Drift length in mm
- `P_bar`: Pressure in bar

# Returns
- Transverse diffusion σ in mm
"""
function electron_transverse_diffusion_Xe_He_mm(L_mm::Real, P_bar::Real)
    return electron_transverse_diffusion_mm(L_mm, P_bar; D_T=D_T_Xe_He)
end

"""
    plot_diffusion_vs_drift(; E::Real=300.0, T::Real=300.0,
                            L_min::Real=10.0, L_max::Real=130.0)

Create a plot of transverse diffusion σ (mm) vs drift length L (mm).

# Arguments
- `E`: Electric field in V/cm (default: 300)
- `T`: Temperature in K (default: 300)
- `L_min`: Minimum drift length in mm (default: 10)
- `L_max`: Maximum drift length in mm (default: 130)

# Returns
- Plot object
"""
function plot_diffusion_vs_drift(; E::Real=300.0, T::Real=300.0,
                                  L_min::Real=10.0, L_max::Real=1300.0)
    L_range = range(L_min, L_max, length=100)
    σ_values = [transverse_diffusion_mm(L, E; T=T) for L in L_range]

    p = plot(L_range, σ_values,
        xlabel = "Drift length L (mm)",
        ylabel = "Transverse diffusion σ (mm)",
        title = "Thermal Ion Diffusion (E = $E V/cm, T = $T K)",
        legend = false,
        linewidth = 1.5,
        linecolor = :black,
        grid = true
    )
    return p
end

"""
    plot_diffusion_vs_drift_multi_E(; E_values::Vector{<:Real}=[200.0, 300.0, 400.0],
                                      T::Real=300.0,
                                      L_min::Real=10.0, L_max::Real=2600.0,
                                      P_bar::Real=15.0,
                                      include_electrons::Bool=true)

Create a plot of transverse diffusion σ (mm) vs drift length L (mm)
for multiple electric field values (ions) and electron diffusion.

# Arguments
- `E_values`: Vector of electric field values in V/cm (default: [200, 300, 400])
- `T`: Temperature in K (default: 300)
- `L_min`: Minimum drift length in mm (default: 10)
- `L_max`: Maximum drift length in mm (default: 2600)
- `P_bar`: Pressure in bar for electron diffusion (default: 15)
- `include_electrons`: Include electron diffusion curves (default: true)

# Returns
- Plot object
"""
function plot_diffusion_vs_drift_multi_E(; E_values::Vector{<:Real}=[200.0, 300.0, 400.0],
                                          T::Real=300.0,
                                          L_min::Real=10.0, L_max::Real=2600.0,
                                          P_bar::Real=15.0,
                                          include_electrons::Bool=true)
    L_range = range(L_min, L_max, length=100)
    colors_ions = [:blue, :red, :green]

    p = plot(
        xlabel = "Drift length L (mm)",
        ylabel = "Transverse diffusion σ (mm)",
        title = "Transverse Diffusion: Ions (T = $T K) & Electrons (P = $P_bar bar)",
        legend = :topleft,
        grid = true,
        size = (900, 600)
    )

    # Plot ion diffusion curves for different E values
    for (i, E) in enumerate(E_values)
        σ_values = [transverse_diffusion_mm(L, E; T=T) for L in L_range]
        plot!(p, L_range, σ_values,
            label = "Ions, E = $(Int(E)) V/cm",
            linewidth = 2,
            linecolor = colors_ions[mod1(i, length(colors_ions))]
        )
    end

    # Plot electron diffusion curves
    if include_electrons
        # Pure Xe electrons
        σ_e_Xe = [electron_transverse_diffusion_pure_Xe_mm(L, P_bar) for L in L_range]
        plot!(p, L_range, σ_e_Xe,
            label = "e⁻ pure Xe (D_T = $D_T_Xe_pure)",
            linewidth = 2,
            linecolor = :black,
            linestyle = :dash
        )

        # Xe/He mixture electrons
        σ_e_XeHe = [electron_transverse_diffusion_Xe_He_mm(L, P_bar) for L in L_range]
        plot!(p, L_range, σ_e_XeHe,
            label = "e⁻ Xe/He 90/10 (D_T = $D_T_Xe_He)",
            linewidth = 2,
            linecolor = :purple,
            linestyle = :dash
        )
    end

    return p
end

"""
    plot_diffusion_anticorrelation(; E::Real=200.0, T::Real=300.0,
                                     P_bar::Real=15.0, L_max::Real=2600.0)

Create a plot showing anti-correlation between ion and electron diffusion.

When ions drift from 0 to L_max, electrons drift from L_max to 0 (opposite directions).
This plot shows both on the same graph with dual x-axes.

# Arguments
- `E`: Electric field in V/cm for ions (default: 200)
- `T`: Temperature in K (default: 300)
- `P_bar`: Pressure in bar (default: 15)
- `L_max`: Maximum drift length in mm (default: 2600)

# Returns
- Plot object
"""
function plot_diffusion_anticorrelation(; E::Real=200.0, T::Real=300.0,
                                          P_bar::Real=15.0, L_max::Real=2600.0)
    n_points = 100
    # Ion drift: 0 → L_max
    L_ion = range(0.0, L_max, length=n_points)
    # Electron drift: L_max → 0 (opposite direction)
    L_electron = range(L_max, 0.0, length=n_points)

    # Compute diffusion values
    # For ions at position L_ion, they have drifted L_ion
    σ_ion = [L > 0 ? transverse_diffusion_mm(L, E; T=T) : 0.0 for L in L_ion]

    # For electrons at the same ion position, electrons have drifted (L_max - L_ion)
    # which equals L_electron
    σ_e_Xe = [L > 0 ? electron_transverse_diffusion_pure_Xe_mm(L, P_bar) : 0.0 for L in L_electron]
    σ_e_XeHe = [L > 0 ? electron_transverse_diffusion_Xe_He_mm(L, P_bar) : 0.0 for L in L_electron]

    # Create plot with ion drift length as x-axis
    p = plot(
        xlabel = "Ion drift length (mm)",
        ylabel = "Transverse diffusion σ (mm)",
        title = "Ion-Electron Diffusion Anti-correlation\n(E = $E V/cm, P = $P_bar bar, T = $T K)",
        legend = :top,
        grid = true,
        size = (900, 600),
        top_margin = 15Plots.mm,
        bottom_margin = 10Plots.mm
    )

    # Plot ion diffusion (increases with ion drift length)
    plot!(p, L_ion, σ_ion,
        label = "Ions (E = $(Int(E)) V/cm)",
        linewidth = 2.5,
        linecolor = :blue
    )

    # Plot electron diffusion (decreases with ion drift length, since electrons drift opposite)
    plot!(p, L_ion, σ_e_Xe,
        label = "e⁻ pure Xe (D_T = $D_T_Xe_pure)",
        linewidth = 2.5,
        linecolor = :black,
        linestyle = :dash
    )

    plot!(p, L_ion, σ_e_XeHe,
        label = "e⁻ Xe/He 90/10 (D_T = $D_T_Xe_He)",
        linewidth = 2.5,
        linecolor = :red,
        linestyle = :dash
    )

    # Add secondary x-axis labels for electron drift length at top
    # We'll add tick annotations at the top
    electron_ticks = [0, 500, 1000, 1500, 2000, 2600]
    ion_positions = [L_max - et for et in electron_ticks]

    # Add electron scale annotation at top
    y_max = maximum([maximum(σ_ion), maximum(σ_e_Xe), maximum(σ_e_XeHe)])
    y_top = y_max * 1.02

    # Add text label for electron axis
    annotate!(p, L_max/2, y_top + 0.15, text("Electron drift length (mm)", 10, :center))

    # Add electron tick labels at top
    for (ion_pos, e_tick) in zip(ion_positions, electron_ticks)
        annotate!(p, ion_pos, y_top + 0.05, text("$e_tick", 8, :center, :gray))
    end

    # Add vertical line at midpoint to show crossing
    vline!(p, [L_max/2], linestyle=:dot, linecolor=:gray, label=nothing, alpha=0.5)

    return p
end

"""
    molar_volume_vdw(P_Pa::Real, T::Real, a::Real, b::Real)

Solve van der Waals equation for molar volume:
    (P + a/V²)(V - b) = RT

This is a cubic equation in V. We solve it numerically using Newton-Raphson.

# Arguments
- `P_Pa`: Pressure in Pa
- `T`: Temperature in K
- `a`: van der Waals constant a (Pa·m⁶/mol²)
- `b`: van der Waals constant b (m³/mol)

# Returns
- Molar volume in m³/mol
"""
function molar_volume_vdw(P_Pa::Real, T::Real, a::Real, b::Real)
    # Initial guess from ideal gas law
    V = R_gas * T / P_Pa

    # Newton-Raphson iteration
    for _ in 1:50
        # f(V) = P - RT/(V-b) + a/V² = 0
        f = P_Pa - R_gas * T / (V - b) + a / V^2
        # f'(V) = RT/(V-b)² - 2a/V³
        df = R_gas * T / (V - b)^2 - 2 * a / V^3

        V_new = V - f / df
        if abs(V_new - V) < 1e-12
            break
        end
        V = V_new
    end

    return V
end

"""
    compressibility_factor(P_bar::Real, T::Real; f_Xe::Real=1.0, f_He::Real=0.0)

Compute the compressibility factor Z = PV/(nRT) for a Xe/He mixture.

Uses van der Waals equation with mixing rules for the constants.

# Arguments
- `P_bar`: Pressure in bar
- `T`: Temperature in K
- `f_Xe`: Fraction of Xenon (default: 1.0)
- `f_He`: Fraction of Helium (default: 0.0)

# Returns
- Compressibility factor Z (dimensionless)
"""
function compressibility_factor(P_bar::Real, T::Real; f_Xe::Real=1.0, f_He::Real=0.0)
    P_Pa = P_bar * 1e5

    # Mixing rules for van der Waals constants (linear mixing for simplicity)
    # For a more accurate treatment, use: a_mix = Σᵢ Σⱼ xᵢxⱼ√(aᵢaⱼ)
    a_mix = f_Xe^2 * a_Xe + f_He^2 * a_He + 2 * f_Xe * f_He * sqrt(a_Xe * a_He)
    b_mix = f_Xe * b_Xe + f_He * b_He

    # Get molar volume from van der Waals
    V_m = molar_volume_vdw(P_Pa, T, a_mix, b_mix)

    # Ideal gas molar volume
    V_ideal = R_gas * T / P_Pa

    # Z = V_real / V_ideal
    return V_m / V_ideal
end

"""
    gas_mass_kg(length_m::Real, diameter_m::Real, P_bar::Real;
                f_Xe::Real=1.0, f_He::Real=0.0, T::Real=300.0, real_gas::Bool=true)

Compute the total mass of gas in a cylindrical TPC.

# Arguments
- `length_m`: TPC length in meters
- `diameter_m`: TPC diameter in meters
- `P_bar`: Pressure in bar
- `f_Xe`: Fraction of Xenon (default: 1.0 for pure Xe)
- `f_He`: Fraction of Helium (default: 0.0)
- `T`: Temperature in K (default: 300 K)
- `real_gas`: Use van der Waals equation (default: true). If false, uses ideal gas law.

# Returns
- Total gas mass in kg
"""
function gas_mass_kg(length_m::Real, diameter_m::Real, P_bar::Real;
                     f_Xe::Real=1.0, f_He::Real=0.0, T::Real=300.0, real_gas::Bool=true)
    @assert f_Xe + f_He ≈ 1.0 "Fractions must sum to 1"
    @assert 0 ≤ f_Xe ≤ 1 && 0 ≤ f_He ≤ 1 "Fractions must be between 0 and 1"

    # Calculate TPC volume (cylinder: V = π r² L)
    radius_m = diameter_m / 2.0
    V_m3 = π * radius_m^2 * length_m

    # Convert pressure from bar to Pa (1 bar = 1e5 Pa)
    P_Pa = P_bar * 1e5

    if real_gas
        # Use van der Waals equation
        # Mixing rules for constants
        a_mix = f_Xe^2 * a_Xe + f_He^2 * a_He + 2 * f_Xe * f_He * sqrt(a_Xe * a_He)
        b_mix = f_Xe * b_Xe + f_He * b_He

        # Get molar volume from van der Waals
        V_molar = molar_volume_vdw(P_Pa, T, a_mix, b_mix)

        # Number of moles = total volume / molar volume
        n_mol = V_m3 / V_molar
    else
        # Ideal gas law: n = PV / (RT) [moles]
        n_mol = (P_Pa * V_m3) / (R_gas * T)
    end

    # Average molar mass of mixture (g/mol)
    M_mix = f_Xe * M_Xe + f_He * M_He

    # Total mass: m = n × M (convert g to kg)
    mass_g = n_mol * M_mix
    mass_kg = mass_g / 1000.0

    return mass_kg
end

"""
    main()

Run all calculations and generate plots.
"""
function main()
    # Parameters from the problem
    f_Xe = 0.90    # 90% Xenon
    f_He = 0.10    # 10% Helium
    P = 15.0       # bar
    E = 300.0      # V/cm
    T = 300.0      # K

    μ_stp = mobility_mixture_stp(f_Xe, f_He)
    μ_P = mobility_at_pressure(μ_stp, P)
    v_d = drift_velocity(f_Xe, f_He, P, E)
    v_d_mm_μs = drift_velocity_mm_μs(f_Xe, f_He, P, E)

    println("\n=== Xe2+ Drift Velocity in Xe/He Mixture ===")
    println("\nInput Parameters:")
    @printf("  Xe fraction: %.0f%%\n", f_Xe * 100)
    @printf("  He fraction: %.0f%%\n", f_He * 100)
    @printf("  Pressure: %.1f bar\n", P)
    @printf("  Electric field: %.1f V/cm\n", E)

    println("\nMobility Values at STP:")
    @printf("  μ(Xe2+ in Xe):  %.2f cm²/V/s\n", μ_Xe2_STP)
    @printf("  μ(Xe2+ in He):  %.2f cm²/V/s\n", μ_He_STP)
    @printf("  μ(NH4+ in Xe):  %.2f cm²/V/s (interpolated from alkaline fit)\n", μ_NH4_STP)
    @printf("  μ(Ba2+ in Xe):  %.2f cm²/V/s\n", μ_Ba2_STP)

    # Alkaline ion fit results
    println("\n=== Alkaline Ion Mobility Fit ===")
    a, b = linear_fit(
        [ALKALINE_IONS.Na.radius, ALKALINE_IONS.K.radius, ALKALINE_IONS.Rb.radius, ALKALINE_IONS.Cs.radius],
        [ALKALINE_IONS.Na.K0, ALKALINE_IONS.K.K0, ALKALINE_IONS.Rb.K0, ALKALINE_IONS.Cs.K0]
    )
    @printf("  Linear fit: K₀ = %.4f + (%.4f) × r\n", a, b)
    @printf("  NH4+ ionic radius: %.2f Å\n", r_NH4)
    @printf("  NH4+ interpolated mobility: %.4f cm²/V/s\n", μ_NH4_STP)

    println("\nResults for Xe2+ in Xe/He mixture:")
    @printf("  Mixture mobility at STP: %.4f cm²/V/s\n", μ_stp)
    @printf("  Mixture mobility at %.1f bar: %.6f cm²/V/s\n", P, μ_P)
    @printf("  Drift velocity: %.4f cm/s\n", v_d)
    @printf("  Drift velocity: %.6f mm/μs\n", v_d_mm_μs)

    # Drift velocities comparison
    println("\n=== Drift Velocities Comparison (mm/s) ===")
    @printf("\nConditions: P = %.1f bar, E = %.1f V/cm\n", P, E)
    println("-" ^ 40)
    @printf("  Xe2+ in pure Xe:     %.2f mm/s\n", drift_velocity_ion_mm_s(μ_Xe2_STP, P, E))
    @printf("  NH4+ in pure Xe:     %.2f mm/s\n", drift_velocity_ion_mm_s(μ_NH4_STP, P, E))
    @printf("  Ba2+ in pure Xe:     %.2f mm/s\n", drift_velocity_ion_mm_s(μ_Ba2_STP, P, E))
    @printf("  Xe2+ in Xe/He mix:   %.2f mm/s\n", v_d * 10.0)

    # Transverse diffusion examples
    println("\n=== Transverse Diffusion (Thermal Ions) ===")
    @printf("\nTemperature: %.1f K, Electric field: %.1f V/cm\n", T, E)
    println("\nDrift length (mm)  |  σ (mm)")
    println("-" ^ 35)
    for L in [10.0, 50.0, 100.0, 130.0]
        σ = transverse_diffusion_mm(L, E; T=T)
        @printf("     %6.1f        |  %.4f\n", L, σ)
    end

    # Generate plots
    # Gas mass calculation
    println("\n=== TPC Gas Mass ===")
    tpc_length = 2.6    # m
    tpc_diameter = 2.6  # m
    tpc_pressure = 15.0 # bar

    # Compressibility factors
    Z_pure_Xe = compressibility_factor(tpc_pressure, T; f_Xe=1.0, f_He=0.0)
    Z_Xe_He = compressibility_factor(tpc_pressure, T; f_Xe=0.9, f_He=0.1)

    # Masses with real gas (van der Waals)
    mass_pure_Xe_real = gas_mass_kg(tpc_length, tpc_diameter, tpc_pressure; f_Xe=1.0, f_He=0.0, T=T, real_gas=true)
    mass_Xe_He_real = gas_mass_kg(tpc_length, tpc_diameter, tpc_pressure; f_Xe=0.9, f_He=0.1, T=T, real_gas=true)

    # Masses with ideal gas (for comparison)
    mass_pure_Xe_ideal = gas_mass_kg(tpc_length, tpc_diameter, tpc_pressure; f_Xe=1.0, f_He=0.0, T=T, real_gas=false)
    mass_Xe_He_ideal = gas_mass_kg(tpc_length, tpc_diameter, tpc_pressure; f_Xe=0.9, f_He=0.1, T=T, real_gas=false)

    @printf("\nTPC dimensions: %.1f m diameter × %.1f m length\n", tpc_diameter, tpc_length)
    @printf("Pressure: %.1f bar, Temperature: %.1f K\n", tpc_pressure, T)
    @printf("TPC volume: %.2f m³\n", π * (tpc_diameter/2)^2 * tpc_length)

    println("\nCompressibility factors (Z = PV/nRT):")
    @printf("  Pure Xe:       Z = %.4f\n", Z_pure_Xe)
    @printf("  Xe/He (90/10): Z = %.4f\n", Z_Xe_He)

    println("\nGas masses:")
    println("-" ^ 55)
    @printf("                      Ideal Gas    Real Gas (vdW)    Δ\n")
    println("-" ^ 55)
    @printf("  Pure Xe:         %8.2f kg    %8.2f kg    %+.1f%%\n",
            mass_pure_Xe_ideal, mass_pure_Xe_real,
            100*(mass_pure_Xe_real - mass_pure_Xe_ideal)/mass_pure_Xe_ideal)
    @printf("  Xe/He (90/10):   %8.2f kg    %8.2f kg    %+.1f%%\n",
            mass_Xe_He_ideal, mass_Xe_He_real,
            100*(mass_Xe_He_real - mass_Xe_He_ideal)/mass_Xe_He_ideal)

    println("\n=== Generating Plots ===")

    # Create output directory
    script_dir = @__DIR__
    output_dir = joinpath(script_dir, "gasProperties")
    mkpath(output_dir)
    println("  Output directory: $output_dir")

    println("  Alkaline ion mobility fit plot...")
    p0 = plot_alkaline_mobility_fit()
    savefig(p0, joinpath(output_dir, "alkaline_mobility_fit.png"))
    println("  -> Saved to gasProperties/alkaline_mobility_fit.png")

    println("  Diffusion vs drift length plot (single E)...")
    p1 = plot_diffusion_vs_drift(; E=E, T=T, L_max=2600.0)
    savefig(p1, joinpath(output_dir, "diffusion_vs_drift.png"))
    println("  -> Saved to gasProperties/diffusion_vs_drift.png")

    println("  Diffusion vs drift length plot (E = 200, 300, 400 V/cm)...")
    p1b = plot_diffusion_vs_drift_multi_E(; E_values=[200.0, 300.0, 400.0], T=T, L_max=2600.0)
    savefig(p1b, joinpath(output_dir, "diffusion_vs_drift_multi_E.png"))
    println("  -> Saved to gasProperties/diffusion_vs_drift_multi_E.png")

    println("  Drift velocities bar plot...")
    p2 = plot_drift_velocities_bar(P, E; f_Xe=f_Xe, f_He=f_He)
    savefig(p2, joinpath(output_dir, "drift_velocities_comparison.png"))
    println("  -> Saved to gasProperties/drift_velocities_comparison.png")

    println("  Drift time vs distance plot...")
    p3 = plot_drift_time_vs_distance(P, E)
    savefig(p3, joinpath(output_dir, "drift_time_vs_distance.png"))
    println("  -> Saved to gasProperties/drift_time_vs_distance.png")

    println("  Arrival time difference plot...")
    p4 = plot_arrival_time_difference(P, E)
    savefig(p4, joinpath(output_dir, "arrival_time_difference.png"))
    println("  -> Saved to gasProperties/arrival_time_difference.png")

    println("  Ion-electron diffusion anti-correlation plot...")
    p5 = plot_diffusion_anticorrelation(; E=200.0, T=T, P_bar=P, L_max=2600.0)
    savefig(p5, joinpath(output_dir, "diffusion_anticorrelation.png"))
    println("  -> Saved to gasProperties/diffusion_anticorrelation.png")

    println("\nDone!")
end

# -------------------------------
# Command line execution
# -------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
