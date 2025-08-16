#!/usr/bin/env julia

"""
HD5t Detector Script

This script creates the HD5t detector configuration for background calculations.
The detector includes copper shielding, PTFE lining, and a central cathode grid.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JHD5t
using Unitful
using Printf
using Plots
using DataFrames
using CSV

"""
    hd5t_detector()

Create HD5t detector for background calculations.
    1. Steel vessel not included (radioactivity from vessel assumed to be shielded by copper)
    2. Copper thickness of 3 cm. This is the effective mass that shoots into detector, due to self-shielding.
    3. PTFE with a thickness of 5 mm in barrel and end-caps.
    4. A central grid of steel with a thickness of 1 mm.
    5. Place all volumes centered in a cylindrical envelope of R=250 cm, L=500 cm
"""
function hd5t_detector()
    
    println("Creating HD5t detector components...")
    println("="^60)
    
    # Define the envelope (mother volume) for all components
    envelope = CylindricalEnvelope(250.0u"cm", 500.0u"cm")
    
    # 1. Barrel Shield (BSL) - Copper shell of Rin = 200 cm, Rout = 203 cm, L=400 cm
    println("1. Barrel Shield (BSL) - Copper")
    bsl_shell = CylinderShell(200.0u"cm", 203.0u"cm", 400.0u"cm")
    bsl = PhysicalCylindricalShell(bsl_shell, copper)
    bsl_placed = PlacedVolume(bsl, "Barrel Shield")
    
    @printf("   Inner radius: %.1f cm\n", ustrip(u"cm", bsl_shell.Rin))
    @printf("   Outer radius: %.1f cm\n", ustrip(u"cm", bsl_shell.Rout))
    @printf("   Length: %.1f cm\n", ustrip(u"cm", bsl_shell.L))
    @printf("   Thickness: %.1f cm\n", ustrip(u"cm", thickness(bsl_shell)))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(bsl)))
    
    # 2. Endcap Shield Left (ESL) - Copper Cylinder of R=203 cm, L=3 cm
    println("2. Endcap Shield Left (ESL) - Copper")
    esl_cyl = Cylinder(203.0u"cm", 3.0u"cm")
    esl = PhysicalCylinder(esl_cyl, copper)
    esl_pos = Position(0.0u"cm", 0.0u"cm", -201.5u"cm")  # Position at left end
    esl_placed = PlacedVolume(esl, esl_pos, "Endcap Shield Left")
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", esl_cyl.R))
    @printf("   Thickness: %.1f cm\n", ustrip(u"cm", esl_cyl.L))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(esl)))
    
    # 3. Endcap Shield Right (ESR) - Copper Cylinder of R=203 cm, L=3 cm
    println("3. Endcap Shield Right (ESR) - Copper")
    esr_cyl = Cylinder(203.0u"cm", 3.0u"cm")
    esr = PhysicalCylinder(esr_cyl, copper)
    esr_pos = Position(0.0u"cm", 0.0u"cm", 201.5u"cm")  # Position at right end
    esr_placed = PlacedVolume(esr, esr_pos, "Endcap Shield Right")
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", esr_cyl.R))
    @printf("   Thickness: %.1f cm\n", ustrip(u"cm", esr_cyl.L))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(esr)))
    
    # 4. Barrel Teflon (BTF) - PTFE shell of Rin=199.5 cm, Rout=200 cm, L=400 cm
    println("4. Barrel Teflon (BTF) - PTFE")
    btf_shell = CylinderShell(199.5u"cm", 200.0u"cm", 400.0u"cm")
    btf = PhysicalCylindricalShell(btf_shell, ptfe)
    btf_placed = PlacedVolume(btf, "Barrel Teflon")
    
    @printf("   Inner radius: %.1f cm\n", ustrip(u"cm", btf_shell.Rin))
    @printf("   Outer radius: %.1f cm\n", ustrip(u"cm", btf_shell.Rout))
    @printf("   Length: %.1f cm\n", ustrip(u"cm", btf_shell.L))
    @printf("   Thickness: %.1f mm\n", ustrip(u"mm", thickness(btf_shell)))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(btf)))
    
    # 5. Endcap Teflon Left (ETL) - PTFE Cylinder of R=199.5 cm, L=5 mm
    println("5. Endcap Teflon Left (ETL) - PTFE")
    etl_cyl = Cylinder(199.5u"cm", 0.5u"cm")  # 5 mm thickness
    etl = PhysicalCylinder(etl_cyl, ptfe)
    etl_pos = Position(0.0u"cm", 0.0u"cm", -200.25u"cm")  # Position at left end (inside copper)
    etl_placed = PlacedVolume(etl, etl_pos, "Endcap Teflon Left")
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", etl_cyl.R))
    @printf("   Thickness: %.1f mm\n", ustrip(u"mm", etl_cyl.L))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(etl)))
    
    # 6. Endcap Teflon Right (ETR) - PTFE Cylinder of R=199.5 cm, L=5 mm
    println("6. Endcap Teflon Right (ETR) - PTFE")
    etr_cyl = Cylinder(199.5u"cm", 0.5u"cm")  # 5 mm thickness
    etr = PhysicalCylinder(etr_cyl, ptfe)
    etr_pos = Position(0.0u"cm", 0.0u"cm", 200.25u"cm")  # Position at right end (inside copper)
    etr_placed = PlacedVolume(etr, etr_pos, "Endcap Teflon Right")
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", etr_cyl.R))
    @printf("   Thickness: %.1f mm\n", ustrip(u"mm", etr_cyl.L))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(etr)))
    
    # 7. Central Cathode (CC) - Fe316Ti Cylinder of R=199.5 cm, L=6.2 um, see mesh.jl

    println("7. Central Cathode (CC) - Fe316Ti")

    cc_cyl = Cylinder(199.5u"cm", 6.2u"μm")  # 6.2 um thickness
    cc = PhysicalCylinder(cc_cyl, fe316ti)
    cc_pos = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")  # Centered in detector
    cc_placed = PlacedVolume(cc, cc_pos, "Central Cathode")
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", cc_cyl.R))
    @printf("   Thickness: %.1f mm\n", ustrip(u"mm", cc_cyl.L))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(cc)))
    
    # 8. Gas Volume - Xenon at 15 bar, 20°C
    println("8. Gas Volume - Xenon")
    gas_cyl = Cylinder(199.5u"cm", 399.9u"cm")  # Inside PTFE lining
    gxe = GXe("rho_1520")  # Create GXe object for xenon at ~15 bar, 20°C
    gas_mat = RadioactiveMaterial(gxe)  # Convert to RadioactiveMaterial with zero radioactivity
    gas = PhysicalCylinder(gas_cyl, gas_mat)
    gas_placed = PlacedVolume(gas, "Gas Volume")
    
    @printf("   Radius: %.1f cm\n", ustrip(u"cm", gas_cyl.R))
    @printf("   Length: %.1f cm\n", ustrip(u"cm", gas_cyl.L))
    @printf("   Volume: %.2f L\n", ustrip(u"L", volume(gas_cyl)))
    @printf("   Mass: %.2f kg\n\n", ustrip(u"kg", mass(gas)))
    
    # Create a custom detector structure with all components
    # Return as a named tuple for easy access
    detector = (
        envelope = envelope,
        bsl = bsl_placed,
        esl = esl_placed,
        esr = esr_placed,
        btf = btf_placed,
        etl = etl_placed,
        etr = etr_placed,
        cc = cc_placed,
        gas = gas_placed,
        # Store all placed volumes for iteration
        volumes = [bsl_placed, esl_placed, esr_placed, btf_placed, 
                  etl_placed, etr_placed, cc_placed, gas_placed]
    )
    
    return detector
end

"""
    analyze_hd5t_detector(detector)

Analyze and display HD5t detector properties.
"""
function analyze_hd5t_detector(detector)
    
    println("\nHD5t Detector Analysis")
    println("="^60)
    
    # Component masses
    println("Component Masses:")
    @printf("   Barrel Shield (Copper):     %.2f kg\n", ustrip(u"kg", mass(detector.bsl)))
    @printf("   Endcap Shield Left (Copper): %.2f kg\n", ustrip(u"kg", mass(detector.esl)))
    @printf("   Endcap Shield Right (Copper): %.2f kg\n", ustrip(u"kg", mass(detector.esr)))
    @printf("   Barrel Teflon (PTFE):       %.2f kg\n", ustrip(u"kg", mass(detector.btf)))
    @printf("   Endcap Teflon Left (PTFE):  %.2f kg\n", ustrip(u"kg", mass(detector.etl)))
    @printf("   Endcap Teflon Right (PTFE): %.2f kg\n", ustrip(u"kg", mass(detector.etr)))
    @printf("   Central Cathode (Fe316Ti):  %.2f kg\n", ustrip(u"kg", mass(detector.cc)))
    @printf("   Gas Volume (Xenon):         %.2f kg\n", ustrip(u"kg", mass(detector.gas)))
    println("   " * "-"^40)
    
    total_mass = sum(mass(vol) for vol in detector.volumes)
    @printf("   Total Mass:                 %.2f kg\n\n", ustrip(u"kg", total_mass))
    
    # Activities for background calculations
    println("Material Activities (Bi-214):")
    @printf("   Copper:   %.1f μBq/kg\n", ustrip(u"μBq/kg", a_bi214(detector.bsl.volume)))
    @printf("   PTFE:     %.1f μBq/kg\n", ustrip(u"μBq/kg", a_bi214(detector.btf.volume)))
    @printf("   Fe316Ti:  %.1f mBq/kg\n", ustrip(u"mBq/kg", a_bi214(detector.cc.volume)))
    @printf("   Xenon:    %.3e Bq/kg\n\n", ustrip(u"Bq/kg", a_bi214(detector.gas.volume)))
    
    println("Material Activities (Tl-208):")
    @printf("   Copper:   %.1f μBq/kg\n", ustrip(u"μBq/kg", a_tl208(detector.bsl.volume)))
    @printf("   PTFE:     %.1f μBq/kg\n", ustrip(u"μBq/kg", a_tl208(detector.btf.volume)))
    @printf("   Fe316Ti:  %.1f mBq/kg\n", ustrip(u"mBq/kg", a_tl208(detector.cc.volume)))
    @printf("   Xenon:    %.3e Bq/kg\n\n", ustrip(u"Bq/kg", a_tl208(detector.gas.volume)))
    
    # Total activities
    println("Total Component Activities:")
    for (name, vol) in [("Barrel Shield", detector.bsl), 
                        ("Endcaps Shield", detector.esl),
                        ("Barrel Teflon", detector.btf),
                        ("Endcaps Teflon", detector.etl),
                        ("Central Cathode", detector.cc)]
        total_bi214 = mass(vol) * a_bi214(vol.volume)
        total_tl208 = mass(vol) * a_tl208(vol.volume)
        @printf("   %s:\n", name)
        @printf("      Bi-214: %.3f mBq total\n", ustrip(u"mBq", total_bi214))
        @printf("      Tl-208: %.3f mBq total\n", ustrip(u"mBq", total_tl208))
    end
    
    # Geometric summary
    println("\nGeometric Summary:")
    println("   Active volume: π × (199.5 cm)² × 399.9 cm = $(round(π * 199.5^2 * 399.9 / 1000, digits=1)) L")
    println("   Copper shield thickness: 3.0 cm")
    println("   PTFE lining thickness: 5.0 mm")
    println("   Central cathode thickness: 1.0 mm")
    bounds = envelope_bounds(detector.envelope)
    println("   Envelope dimensions: R=$(ustrip(u"cm", detector.envelope.shape.R)) cm, L=$(ustrip(u"cm", detector.envelope.shape.L)) cm")
end

"""
    create_detector_summary_df(detector)

Create a summary DataFrame with detector component masses and activities.

# Returns
- DataFrame with columns: Component, Mass_kg, Bi214_Activity_mBq_kg, Tl208_Activity_mBq_kg, 
  Total_Bi214_mBq, Total_Tl208_mBq
"""
function create_detector_summary_df(detector)
    
    # Define component data
    components = [
        ("Barrel_Shield_Cu", detector.bsl),
        ("Endcap_Shield_Left_Cu", detector.esl),
        ("Endcap_Shield_Right_Cu", detector.esr),
        ("Barrel_Teflon_PTFE", detector.btf),
        ("Endcap_Teflon_Left_PTFE", detector.etl),
        ("Endcap_Teflon_Right_PTFE", detector.etr),
        ("Central_Cathode_Fe316Ti", detector.cc),
        ("Gas_Volume_Xe", detector.gas)
    ]
    
    # Initialize arrays for DataFrame columns
    component_names = String[]
    masses_kg = Float64[]
    bi214_activities_mBq_kg = Float64[]
    tl208_activities_mBq_kg = Float64[]
    total_bi214_mBq = Float64[]
    total_tl208_mBq = Float64[]
    
    # Populate data
    for (name, vol) in components
        push!(component_names, name)
        
        # Mass in kg
        mass_kg = ustrip(u"kg", mass(vol))
        push!(masses_kg, mass_kg)
        
        # Specific activities in mBq/kg
        bi214_act = ustrip(u"mBq/kg", a_bi214(vol.volume))
        tl208_act = ustrip(u"mBq/kg", a_tl208(vol.volume))
        push!(bi214_activities_mBq_kg, bi214_act)
        push!(tl208_activities_mBq_kg, tl208_act)
        
        # Total activities in mBq
        total_bi214 = ustrip(u"mBq", mass(vol) * a_bi214(vol.volume))
        total_tl208 = ustrip(u"mBq", mass(vol) * a_tl208(vol.volume))
        push!(total_bi214_mBq, total_bi214)
        push!(total_tl208_mBq, total_tl208)
    end
    
    # Create DataFrame
    df = DataFrame(
        Component = component_names,
        Mass_kg = masses_kg,
        Bi214_Activity_mBq_kg = bi214_activities_mBq_kg,
        Tl208_Activity_mBq_kg = tl208_activities_mBq_kg,
        Total_Bi214_mBq = total_bi214_mBq,
        Total_Tl208_mBq = total_tl208_mBq
    )
    
    # Add totals row
    totals_row = DataFrame(
        Component = ["TOTAL"],
        Mass_kg = [sum(masses_kg)],
        Bi214_Activity_mBq_kg = [NaN],  # Not meaningful for total
        Tl208_Activity_mBq_kg = [NaN],  # Not meaningful for total
        Total_Bi214_mBq = [sum(total_bi214_mBq)],
        Total_Tl208_mBq = [sum(total_tl208_mBq)]
    )
    
    # Combine main dataframe with totals
    df = vcat(df, totals_row)
    
    return df
end

"""
    format_detector_summary_df(df; round_digits=2)

Format the detector summary DataFrame with rounded values for better display.
"""
function format_detector_summary_df(df; round_digits=2)
    df_formatted = copy(df)
    
    # Round numerical columns
    df_formatted.Mass_kg = round.(df_formatted.Mass_kg, digits=round_digits)
    df_formatted.Total_Bi214_mBq = round.(df_formatted.Total_Bi214_mBq, digits=round_digits)
    df_formatted.Total_Tl208_mBq = round.(df_formatted.Total_Tl208_mBq, digits=round_digits)
    
    # Format activity columns with appropriate precision
    for i in 1:nrow(df_formatted)
        if !isnan(df_formatted.Bi214_Activity_mBq_kg[i])
            if df_formatted.Bi214_Activity_mBq_kg[i] < 0.1
                df_formatted.Bi214_Activity_mBq_kg[i] = round(df_formatted.Bi214_Activity_mBq_kg[i], digits=4)
            elseif df_formatted.Bi214_Activity_mBq_kg[i] < 1.0
                df_formatted.Bi214_Activity_mBq_kg[i] = round(df_formatted.Bi214_Activity_mBq_kg[i], digits=3)
            else
                df_formatted.Bi214_Activity_mBq_kg[i] = round(df_formatted.Bi214_Activity_mBq_kg[i], digits=2)
            end
        end
        
        if !isnan(df_formatted.Tl208_Activity_mBq_kg[i])
            if df_formatted.Tl208_Activity_mBq_kg[i] < 0.1
                df_formatted.Tl208_Activity_mBq_kg[i] = round(df_formatted.Tl208_Activity_mBq_kg[i], digits=4)
            elseif df_formatted.Tl208_Activity_mBq_kg[i] < 1.0
                df_formatted.Tl208_Activity_mBq_kg[i] = round(df_formatted.Tl208_Activity_mBq_kg[i], digits=3)
            else
                df_formatted.Tl208_Activity_mBq_kg[i] = round(df_formatted.Tl208_Activity_mBq_kg[i], digits=2)
            end
        end
    end
    
    return df_formatted
end

"""
    save_detector_summary(detector; filename="hd5t_detector_summary.csv", output_dir=".")

Save the detector summary DataFrame to a CSV file.

# Arguments
- `detector`: The detector object
- `filename`: Name of the CSV file to save
- `output_dir`: Directory to save the file
"""
function save_detector_summary(detector; filename="hd5t_detector_summary.csv", output_dir=".")
    df = create_detector_summary_df(detector)
    
    # Ensure output directory exists
    if !isdir(output_dir)
        mkpath(output_dir)
    end
    
    # Save to CSV
    filepath = joinpath(output_dir, filename)
    CSV.write(filepath, df)
    
    println("Detector summary saved to: $filepath")
    return df
end




"""
    plot_hd5t_projections(detector; save_plots=false, output_dir=".")

Draw XY and XZ projections of the HD5t detector showing different materials with distinct colors.
"""
function plot_hd5t_projections(detector; save_plots=false, output_dir=".")
    
    # Define material colors
    colors = Dict(
        "copper" => :orange,     # Copper shielding
        "ptfe" => :purple,       # PTFE lining
        "steel" => :gray,        # Central cathode
        "xenon" => :lightblue    # Gas volume
    )
    
    # Create XY projection (top view)
    xy_plot = plot(
        title="HD5t Detector - XY Projection (Top View)",
        xlabel="X (cm)", ylabel="Y (cm)",
        aspect_ratio=:equal,
        legend=:topright,
        size=(800, 800)
    )
    
    # Plot concentric circles
    θ = range(0, 2π, length=100)
    
    # Copper shield (outer)
    r_cu_out = 203.0
    r_cu_in = 200.0
    x_cu_out = r_cu_out .* cos.(θ)
    y_cu_out = r_cu_out .* sin.(θ)
    x_cu_in = r_cu_in .* cos.(θ)
    y_cu_in = r_cu_in .* sin.(θ)
    
    plot!(xy_plot, x_cu_out, y_cu_out, linewidth=3, color=colors["copper"], 
          label="Copper Shield", fill=false)
    plot!(xy_plot, x_cu_in, y_cu_in, linewidth=2, color=colors["copper"], 
          label="", fill=false, alpha=0.7)
    
    # PTFE lining
    r_ptfe_out = 200.0
    r_ptfe_in = 199.5
    x_ptfe_out = r_ptfe_out .* cos.(θ)
    y_ptfe_out = r_ptfe_out .* sin.(θ)
    x_ptfe_in = r_ptfe_in .* cos.(θ)
    y_ptfe_in = r_ptfe_in .* sin.(θ)
    
    plot!(xy_plot, x_ptfe_out, y_ptfe_out, linewidth=2, color=colors["ptfe"], 
          label="PTFE Lining", fill=false)
    plot!(xy_plot, x_ptfe_in, y_ptfe_in, linewidth=1, color=colors["ptfe"], 
          label="", fill=false, alpha=0.7)
    
    # Gas volume
    r_gas = 199.5
    x_gas = r_gas .* cos.(θ)
    y_gas = r_gas .* sin.(θ)
    plot!(xy_plot, x_gas, y_gas, linewidth=1, color=colors["xenon"], 
          label="Gas Volume", fill=true, alpha=0.3)
    
    # Central cathode (visible as a line in XY view at z=0)
    plot!(xy_plot, [-r_gas, r_gas], [0, 0], linewidth=2, color=colors["steel"], 
          label="Central Cathode")
    plot!(xy_plot, [0, 0], [-r_gas, r_gas], linewidth=2, color=colors["steel"], 
          label="")
    
    # Add center point
    scatter!(xy_plot, [0], [0], color=:red, markersize=3, label="Center")
    
    # Create XZ projection (side view)
    xz_plot = plot(
        title="HD5t Detector - XZ Projection (Side View)",
        xlabel="X (cm)", ylabel="Z (cm)",
        aspect_ratio=:equal,
        legend=:topright,
        size=(1000, 600)
    )
    
    # Z positions
    z_center = 0.0
    cu_z_min = -203.0
    cu_z_max = 203.0
    ptfe_z_min = -200.0
    ptfe_z_max = 200.0
    gas_z_min = -199.95
    gas_z_max = 199.95
    
    # Copper shield boundaries
    plot!(xz_plot, [-203, 203], [cu_z_min, cu_z_min], linewidth=3, 
          color=colors["copper"], label="Copper Shield")
    plot!(xz_plot, [-203, 203], [cu_z_max, cu_z_max], linewidth=3, 
          color=colors["copper"], label="")
    plot!(xz_plot, [203, 203], [cu_z_min, cu_z_max], linewidth=3, 
          color=colors["copper"], label="")
    plot!(xz_plot, [-203, -203], [cu_z_min, cu_z_max], linewidth=3, 
          color=colors["copper"], label="")
    
    # Inner copper boundary
    plot!(xz_plot, [-200, 200], [ptfe_z_min, ptfe_z_min], linewidth=2, 
          color=colors["copper"], label="", alpha=0.7)
    plot!(xz_plot, [-200, 200], [ptfe_z_max, ptfe_z_max], linewidth=2, 
          color=colors["copper"], label="", alpha=0.7)
    plot!(xz_plot, [200, 200], [ptfe_z_min, ptfe_z_max], linewidth=2, 
          color=colors["copper"], label="", alpha=0.7)
    plot!(xz_plot, [-200, -200], [ptfe_z_min, ptfe_z_max], linewidth=2, 
          color=colors["copper"], label="", alpha=0.7)
    
    # PTFE lining
    plot!(xz_plot, [-200, 200], [ptfe_z_min, ptfe_z_min], linewidth=2, 
          color=colors["ptfe"], label="PTFE Lining")
    plot!(xz_plot, [-200, 200], [ptfe_z_max, ptfe_z_max], linewidth=2, 
          color=colors["ptfe"], label="")
    plot!(xz_plot, [-199.5, 199.5], [gas_z_min, gas_z_min], linewidth=1, 
          color=colors["ptfe"], label="", alpha=0.7)
    plot!(xz_plot, [-199.5, 199.5], [gas_z_max, gas_z_max], linewidth=1, 
          color=colors["ptfe"], label="", alpha=0.7)
    
    # PTFE endcaps
    plot!(xz_plot, [-199.5, 199.5], [-200.25, -200.25], linewidth=4, 
          color=colors["ptfe"], label="")
    plot!(xz_plot, [-199.5, 199.5], [200.25, 200.25], linewidth=4, 
          color=colors["ptfe"], label="")
    
    # Gas volume (filled)
    plot!(xz_plot, [-199.5, 199.5, 199.5, -199.5, -199.5], 
          [gas_z_min, gas_z_min, gas_z_max, gas_z_max, gas_z_min],
          linewidth=1, color=colors["xenon"], label="Gas Volume", 
          fill=true, alpha=0.3)
    
    # Central cathode
    plot!(xz_plot, [-199.5, 199.5], [0, 0], linewidth=3, 
          color=colors["steel"], label="Central Cathode")
    
    # Add center line
    plot!(xz_plot, [0, 0], [cu_z_min-10, cu_z_max+10], 
          color=:red, linestyle=:dash, alpha=0.5, label="Detector Axis")
    
    # Save plots if requested
    if save_plots
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        savefig(xy_plot, joinpath(output_dir, "hd5t_xy_projection.png"))
        savefig(xz_plot, joinpath(output_dir, "hd5t_xz_projection.png"))
        println("Plots saved to: $output_dir")
    end
    
    return (xy_plot, xz_plot)
end

"""
    plot_hd5t_combined(detector; save_plot=false, output_dir=".")

Create a combined view of the HD5t detector projections.
"""
function plot_hd5t_combined(detector; save_plot=false, output_dir=".")
    xy_plot, xz_plot = plot_hd5t_projections(detector, save_plots=false)
    
    combined_plot = plot(xy_plot, xz_plot, 
                        layout=(1,2), 
                        size=(1600, 800),
                        plot_title="HD5t Detector - Multiple Projections")
    
    if save_plot
        if !isdir(output_dir)
            mkpath(output_dir)
        end
        savefig(combined_plot, joinpath(output_dir, "hd5t_combined.png"))
        println("Combined plot saved to: $output_dir")
    end
    
    return combined_plot
end

function main()
    """Main function to run the HD5t detector analysis"""
    
    println("HD5t Detector Configuration")
    println("="^60)
    println("Creating detector for background calculations...")
    println()
    
    # Create the detector
    detector = hd5t_detector()
    
    # Analyze the detector
    analyze_hd5t_detector(detector)
    
    # Create detector visualizations
    println("\nCreating detector visualizations...")
    println("="^60)
    
    # Generate XY and XZ projections
    xy_plot, xz_plot = plot_hd5t_projections(detector, save_plots=true, output_dir="hd5t_plots")
    
    # Generate combined plot
    combined_plot = plot_hd5t_combined(detector, save_plot=true, output_dir="hd5t_plots")
    
    println("✅ HD5t detector plots created and saved to: hd5t_plots/")
    println("   - hd5t_xy_projection.png")
    println("   - hd5t_xz_projection.png") 
    println("   - hd5t_combined.png")
    
    # Create and save summary DataFrame
    println("\nCreating detector summary DataFrame...")
    println("="^60)
    df_summary = save_detector_summary(detector, output_dir=".")
    
    # Display the formatted summary table
    println("\nDetector Summary Table:")
    df_formatted = format_detector_summary_df(df_summary)
    println(df_formatted)
    
    println("\nHD5t detector analysis completed successfully! 🎉")
end

# Run the script if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end