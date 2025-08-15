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
    @printf("   Density: %.2f g/cm³\n", ustrip(u"g/cm^3", vessel.bst.material.ρ))
    @printf("   Attenuation length: %.1f cm\n", ustrip(u"cm", att_length_bst(vessel)))
    @printf("   Bi-214 activity: %.1f mBq/kg\n", ustrip(u"mBq/kg", a_bi214_bst(vessel)))
    @printf("   Tl-208 activity: %.1f mBq/kg\n\n", ustrip(u"mBq/kg", a_tl208_bst(vessel)))
    
    # Material properties (shield)
    println("Barrel Shield Material Properties:")
    @printf("   Density: %.2f g/cm³\n", ustrip(u"g/cm^3", vessel.bsl.material.ρ))
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
    
    println("Demo completed successfully! 🎉")
end

# Run the script if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end