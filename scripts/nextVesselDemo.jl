#!/usr/bin/env julia

"""
JHD5t NextVessel Demo Script

This script demonstrates the creation and analysis of a NextVessel object,
showing the dimensions and masses of all its components.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JHD5t
using Unitful
using Printf
using Plots

function create_hd5t_vessel()
    """Create HD5t vessel"""
    
    println("Creating HD5t components...")
    println("="^60)
    
    # Barrel Structure (BST) - Steel shell
    println("1. Barrel Structure (BST) - Steel")
    bst_shell = CylinderShell(215.0u"cm", 220.0u"cm", 450.0u"cm")
    bst_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 1.0u"mBq/kg")
    bst = PhysicalCylindricalShell(bst_shell, bst_mat)
    
    @printf("   Inner radius: %.1f cm\n", ustrip(u"cm", bst_shell.Rin))
    @printf("   Outer radius: %.1f cm\n", ustrip(u"cm", bst_shell.Rout))
    @printf("   Length: %.1f cm\n", ustrip(u"cm", bst_shell.L))
    @printf("   Thickness: %.1f cm\n", ustrip(u"cm", thickness(bst_shell)))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(bst)))
    
    # Barrel Shield (BSL) - Copper shell
    println("2. Barrel Shield (BSL) - Copper")
    bsl_shell = CylinderShell(200.0u"cm", 215.0u"cm", 400.0u"cm")
    bsl_mat = RadioactiveMaterial(8.96u"g/cm^3", 1.5u"cm", 10.0u"μBq/kg", 10.0u"μBq/kg")
    bsl = PhysicalCylindricalShell(bsl_shell, bsl_mat)
    
    @printf("   Inner radius: %.1f cm\n", ustrip(u"cm", bsl_shell.Rin))
    @printf("   Outer radius: %.1f cm\n", ustrip(u"cm", bsl_shell.Rout))
    @printf("   Length: %.1f cm\n", ustrip(u"cm", bsl_shell.L))
    @printf("   Thickness: %.1f cm\n", ustrip(u"cm", thickness(bsl_shell)))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(bsl)))
    
    # Endcap Structure (EST) - Steel endcaps
    println("3. Endcap Structure (EST) - Steel")
    est_cyl = Cylinder(220.0u"cm", 5.0u"cm")
    est_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 1.0u"mBq/kg")
    est = PhysicalCylinder(est_cyl, est_mat)
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", est_cyl.R))
    @printf("   Thickness: %.1f cm\n", ustrip(u"cm", est_cyl.L))
    @printf("   Volume: %.2f L\n", ustrip(u"L", volume(est_cyl)))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(est)))
    
    # Endcap Shield (ESL) - Copper endcaps
    println("4. Endcap Shield (ESL) - Copper")
    esl_cyl = Cylinder(215.0u"cm", 15.0u"cm")
    esl_mat = RadioactiveMaterial(8.96u"g/cm^3", 1.5u"cm", 10.0u"μBq/kg", 10.0u"μBq/kg")
    esl = PhysicalCylinder(esl_cyl, esl_mat)
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", esl_cyl.R))
    @printf("   Thickness: %.1f cm\n", ustrip(u"cm", esl_cyl.L))
    @printf("   Volume: %.2f L\n", ustrip(u"L", volume(esl_cyl)))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(esl)))
    
    # Gas Volume - Xenon at 10 bar, 20°C
    println("5. Gas Volume - Xenon")
    gas_cyl = Cylinder(200.0u"cm", 400.0u"cm")
    gxe = GXe("rho_1520")  # Create GXe object for xenon at ~15 bar, 20°C
    gas_mat = RadioactiveMaterial(gxe)  # Convert to RadioactiveMaterial with zero radioactivity
    gas = PhysicalCylinder(gas_cyl, gas_mat)
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", gas_cyl.R))
    @printf("   Length: %.1f cm\n", ustrip(u"cm", gas_cyl.L))
    @printf("   Volume: %.2f L\n", ustrip(u"L", volume(gas_cyl)))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(gas)))
    
    # Create NextVessel
    vessel = NextVessel(bst, bsl, est, esl, gas)
    
    return vessel
end

function analyze_next_vessel(vessel::NextVessel)
    """Analyze and display NextVessel properties"""
    
    println("NextVessel Analysis")
    println("="^60)
    
    # Component masses
    println("Component Masses:")
    @printf("   Barrel Structure (BST): %.2f kg\n", ustrip(u"kg", mass_bst(vessel)))
    @printf("   Barrel Shield (BSL):    %.2f kg\n", ustrip(u"kg", mass_bsl(vessel)))
    @printf("   Endcap Structure (EST): %.2f kg\n", ustrip(u"kg", mass_est(vessel)))
    @printf("   Endcap Shield (ESL):    %.2f kg\n", ustrip(u"kg", mass_esl(vessel)))
    @printf("   Gas Volume:             %.2f kg\n", ustrip(u"kg", mass_gas(vessel)))
    println("   " * "-"^40)
    @printf("   Total Mass:             %.2f kg\n\n", ustrip(u"kg", mass(vessel)))
    
    # Material properties (structure)
    println("Barrel Structure Material Properties:")
    @printf("   Density: %.2f g/cm³\n", ustrip(u"g/cm^3", vessel.bst.volume.material.ρ))
    @printf("   Attenuation length: %.1f cm\n", ustrip(u"cm", att_length_bst(vessel)))
    @printf("   Bi-214 activity: %.1f mBq/kg\n", ustrip(u"mBq/kg", a_bi214_bst(vessel)))
    @printf("   Tl-208 activity: %.1f mBq/kg\n\n", ustrip(u"mBq/kg", a_tl208_bst(vessel)))
    
    # Material properties (shield)
    println("Barrel Shield Material Properties:")
    @printf("   Density: %.2f g/cm³\n", ustrip(u"g/cm^3", vessel.bsl.volume.material.ρ))
    @printf("   Attenuation length: %.1f cm\n", ustrip(u"cm", att_length_bsl(vessel)))
    @printf("   Bi-214 activity: %.1f μBq/kg\n", ustrip(u"μBq/kg", a_bi214_bsl(vessel)))
    @printf("   Tl-208 activity: %.1f μBq/kg\n\n", ustrip(u"μBq/kg", a_tl208_bsl(vessel)))
    
    # Geometric summary
    println("Geometric Summary:")
    println("   Inner active volume: π × (100 cm)² × 200 cm = $(round(π * 100^2 * 200 / 1000, digits=1)) L")
    println("   Barrel wall thickness: 5.0 cm (2 cm steel + 3 cm copper)")
    println("   Endcap thickness: 5.0 cm (2 cm steel + 3 cm copper)")
    println("   Overall dimensions: Ø210 cm × 210 cm")
end

"""
    plot_detector_projections(vessel::NextVessel; save_plots=false, output_dir=".")

Draw XY and XZ projections of the HD5t detector showing different materials with distinct colors.

# Arguments
- `vessel::NextVessel`: The detector vessel to plot
- `save_plots::Bool=false`: Whether to save plots to disk
- `output_dir::String="."`: Directory to save plots (if save_plots=true)

# Returns
- Tuple of (xy_plot, xz_plot): The two projection plots
"""
function plot_detector_projections(vessel::NextVessel; save_plots=false, output_dir=".")
    
    # Extract dimensions (convert to cm for plotting)
    bst_rin = ustrip(u"cm", vessel.bst.volume.shell.Rin)   # 215 cm
    bst_rout = ustrip(u"cm", vessel.bst.volume.shell.Rout) # 220 cm
    bst_length = ustrip(u"cm", vessel.bst.volume.shell.L)  # 450 cm
    
    bsl_rin = ustrip(u"cm", vessel.bsl.volume.shell.Rin)   # 200 cm
    bsl_rout = ustrip(u"cm", vessel.bsl.volume.shell.Rout) # 215 cm
    bsl_length = ustrip(u"cm", vessel.bsl.volume.shell.L)  # 400 cm
    
    gas_r = ustrip(u"cm", vessel.gas.volume.shell.R)       # 200 cm
    gas_length = ustrip(u"cm", vessel.gas.volume.shell.L)  # 400 cm
    
    est_r = ustrip(u"cm", vessel.est.volume.shell.R)       # 220 cm
    est_thickness = ustrip(u"cm", vessel.est.volume.shell.L) # 5 cm
    
    esl_r = ustrip(u"cm", vessel.esl.volume.shell.R)       # 215 cm
    esl_thickness = ustrip(u"cm", vessel.esl.volume.shell.L) # 15 cm
    
    # Define material colors
    colors = Dict(
        "steel" => :gray,      # Barrel Structure (BST) + Endcap Structure (EST)
        "copper" => :orange,   # Barrel Shield (BSL) + Endcap Shield (ESL)
        "xenon" => :lightblue  # Gas volume
    )
    
    # Create XY projection (top view)
    xy_plot = plot(
        title="NextVessel - XY Projection (Top View)",
        xlabel="X (cm)", ylabel="Y (cm)",
        aspect_ratio=:equal,
        legend=:topright,
        size=(800, 800)
    )
    
    # Plot concentric circles for barrel components (XY view)
    θ = range(0, 2π, length=100)
    
    # Barrel Structure (Steel) - outermost
    x_bst = bst_rout .* cos.(θ)
    y_bst = bst_rout .* sin.(θ)
    plot!(xy_plot, x_bst, y_bst, linewidth=3, color=colors["steel"], 
          label="Barrel Structure (Steel)", fill=false)
    
    # Inner edge of Barrel Structure
    x_bst_in = bst_rin .* cos.(θ)
    y_bst_in = bst_rin .* sin.(θ)
    plot!(xy_plot, x_bst_in, y_bst_in, linewidth=2, color=colors["steel"], 
          label="", fill=false, alpha=0.7)
    
    # Barrel Shield (Copper)
    x_bsl = bsl_rout .* cos.(θ)
    y_bsl = bsl_rout .* sin.(θ)
    plot!(xy_plot, x_bsl, y_bsl, linewidth=3, color=colors["copper"], 
          label="Barrel Shield (Copper)", fill=false)
    
    # Inner edge of Barrel Shield
    x_bsl_in = bsl_rin .* cos.(θ)
    y_bsl_in = bsl_rin .* sin.(θ)
    plot!(xy_plot, x_bsl_in, y_bsl_in, linewidth=2, color=colors["copper"], 
          label="", fill=false, alpha=0.7)
    
    # Gas volume (Xenon)
    x_gas = gas_r .* cos.(θ)
    y_gas = gas_r .* sin.(θ)
    plot!(xy_plot, x_gas, y_gas, linewidth=2, color=colors["xenon"], 
          label="Gas Volume (Xenon)", fill=true, alpha=0.3)
    
    # Add center point
    scatter!(xy_plot, [0], [0], color=:red, markersize=3, label="Center")
    
    # Create XZ projection (side view)
    xz_plot = plot(
        title="NextVessel - XZ Projection (Side View)",
        xlabel="X (cm)", ylabel="Z (cm)",
        aspect_ratio=:equal,
        legend=:topright,
        size=(1000, 600)
    )
    
    # Calculate Z positions (centered at origin)
    z_center = 0.0
    
    # Barrel Structure extends ±bst_length/2
    bst_z_min = z_center - bst_length/2
    bst_z_max = z_center + bst_length/2
    
    # Barrel Shield extends ±bsl_length/2  
    bsl_z_min = z_center - bsl_length/2
    bsl_z_max = z_center + bsl_length/2
    
    # Gas volume extends ±gas_length/2
    gas_z_min = z_center - gas_length/2
    gas_z_max = z_center + gas_length/2
    
    # Endcap positions
    est_z_left = bst_z_min - est_thickness
    est_z_right = bst_z_max + est_thickness
    esl_z_left = est_z_left + est_thickness
    esl_z_right = est_z_right - est_thickness
    
    # Plot barrel components (XZ view)
    # Barrel Structure (Steel) - top and bottom edges
    plot!(xz_plot, [-bst_rout, bst_rout], [bst_z_min, bst_z_min], 
          linewidth=3, color=colors["steel"], label="Barrel Structure (Steel)")
    plot!(xz_plot, [-bst_rout, bst_rout], [bst_z_max, bst_z_max], 
          linewidth=3, color=colors["steel"], label="")
    plot!(xz_plot, [-bst_rin, bst_rin], [bst_z_min, bst_z_min], 
          linewidth=2, color=colors["steel"], label="", alpha=0.7)
    plot!(xz_plot, [-bst_rin, bst_rin], [bst_z_max, bst_z_max], 
          linewidth=2, color=colors["steel"], label="", alpha=0.7)
    
    # Barrel Structure sides
    plot!(xz_plot, [bst_rout, bst_rout], [bst_z_min, bst_z_max], 
          linewidth=3, color=colors["steel"], label="")
    plot!(xz_plot, [-bst_rout, -bst_rout], [bst_z_min, bst_z_max], 
          linewidth=3, color=colors["steel"], label="")
    plot!(xz_plot, [bst_rin, bst_rin], [bst_z_min, bst_z_max], 
          linewidth=2, color=colors["steel"], label="", alpha=0.7)
    plot!(xz_plot, [-bst_rin, -bst_rin], [bst_z_min, bst_z_max], 
          linewidth=2, color=colors["steel"], label="", alpha=0.7)
    
    # Barrel Shield (Copper)
    plot!(xz_plot, [-bsl_rout, bsl_rout], [bsl_z_min, bsl_z_min], 
          linewidth=3, color=colors["copper"], label="Barrel Shield (Copper)")
    plot!(xz_plot, [-bsl_rout, bsl_rout], [bsl_z_max, bsl_z_max], 
          linewidth=3, color=colors["copper"], label="")
    plot!(xz_plot, [-bsl_rin, bsl_rin], [bsl_z_min, bsl_z_min], 
          linewidth=2, color=colors["copper"], label="", alpha=0.7)
    plot!(xz_plot, [-bsl_rin, bsl_rin], [bsl_z_max, bsl_z_max], 
          linewidth=2, color=colors["copper"], label="", alpha=0.7)
    
    # Barrel Shield sides
    plot!(xz_plot, [bsl_rout, bsl_rout], [bsl_z_min, bsl_z_max], 
          linewidth=3, color=colors["copper"], label="")
    plot!(xz_plot, [-bsl_rout, -bsl_rout], [bsl_z_min, bsl_z_max], 
          linewidth=3, color=colors["copper"], label="")
    plot!(xz_plot, [bsl_rin, bsl_rin], [bsl_z_min, bsl_z_max], 
          linewidth=2, color=colors["copper"], label="", alpha=0.7)
    plot!(xz_plot, [-bsl_rin, -bsl_rin], [bsl_z_min, bsl_z_max], 
          linewidth=2, color=colors["copper"], label="", alpha=0.7)
    
    # Gas volume (Xenon) - filled rectangle
    plot!(xz_plot, [-gas_r, gas_r, gas_r, -gas_r, -gas_r], 
          [gas_z_min, gas_z_min, gas_z_max, gas_z_max, gas_z_min],
          linewidth=2, color=colors["xenon"], label="Gas Volume (Xenon)", 
          fill=true, alpha=0.3)
    
    # Endcap Structure (Steel)
    plot!(xz_plot, [-est_r, est_r], [est_z_left, est_z_left], 
          linewidth=4, color=colors["steel"], label="Endcap Structure (Steel)")
    plot!(xz_plot, [-est_r, est_r], [est_z_right, est_z_right], 
          linewidth=4, color=colors["steel"], label="")
    
    # Endcap Shield (Copper)
    plot!(xz_plot, [-esl_r, esl_r], [esl_z_left, esl_z_left], 
          linewidth=4, color=colors["copper"], label="Endcap Shield (Copper)")
    plot!(xz_plot, [-esl_r, esl_r], [esl_z_right, esl_z_right], 
          linewidth=4, color=colors["copper"], label="")
    
    # Add center line
    plot!(xz_plot, [0, 0], [est_z_left-10, est_z_right+10], 
          color=:red, linestyle=:dash, alpha=0.5, label="Detector Axis")
    
    # Save plots if requested
    if save_plots
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        savefig(xy_plot, joinpath(output_dir, "next_vessel_xy_projection.png"))
        savefig(xz_plot, joinpath(output_dir, "next_vessel_xz_projection.png"))
        println("Plots saved to: $output_dir")
    end
    
    return (xy_plot, xz_plot)
end

"""
    plot_detector_3d_schematic(vessel::NextVessel; save_plot=false, output_dir=".")

Create a simplified 3D schematic representation of the detector.
"""
function plot_detector_3d_schematic(vessel::NextVessel; save_plot=false, output_dir=".")
    # This would require a 3D plotting backend like PlotlyJS
    # For now, return the 2D projections arranged side by side
    xy_plot, xz_plot = plot_detector_projections(vessel, save_plots=false)
    
    combined_plot = plot(xy_plot, xz_plot, 
                        layout=(1,2), 
                        size=(1600, 800),
                        plot_title="NextVessel - Multiple Projections")
    
    if save_plot
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        savefig(combined_plot, joinpath(output_dir, "next_vessel_combined.png"))
        println("Combined plot saved to: $output_dir")
    end
    
    return combined_plot
end

function main()
    """Main function to run the demo"""
    
    println("JHD5t NextVessel Demonstration")
    println("="^60)
    println("Creating and analyzing a Next-like detector vessel...")
    println()
    
    # Create the vessel
    vessel = create_hd5t_vessel()
    
    # Analyze the vessel
    analyze_next_vessel(vessel)
    
    # Create detector visualizations
    println("\nCreating detector visualizations...")
    println("="^60)
    
    # Generate XY and XZ projections
    xy_plot, xz_plot = plot_detector_projections(vessel, save_plots=true, output_dir="detector_plots")
    
    # Generate combined plot
    combined_plot = plot_detector_3d_schematic(vessel, save_plot=true, output_dir="detector_plots")
    
    println("✅ Detector plots created and saved to: detector_plots/")
    println("   - next_vessel_xy_projection.png")
    println("   - next_vessel_xz_projection.png") 
    println("   - next_vessel_combined.png")
    
    println("\nDemo completed successfully! 🎉")
end

# Run the script if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end