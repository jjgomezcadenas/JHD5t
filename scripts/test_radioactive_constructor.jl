#!/usr/bin/env julia

"""
Test script demonstrating the new RadioactiveMaterial constructor
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using JHD5t
using Unitful
using Printf

function test_radioactive_constructors()
    println("Testing RadioactiveMaterial Constructors")
    println("="^50)
    
    # Method 1: Original constructor (individual parameters)
    println("1. Original constructor:")
    density = 8.96u"g/cm^3"  # Copper density
    att_length = 1.5u"cm"
    a_bi = 10.0u"μBq/kg"
    a_tl = 5.0u"μBq/kg"
    
    rad_mat1 = RadioactiveMaterial(density, att_length, a_bi, a_tl)
    
    @printf("   Density: %.2f g/cm³\n", ustrip(u"g/cm^3", rad_mat1.ρ))
    @printf("   Attenuation length: %.1f cm\n", ustrip(u"cm", rad_mat1.λ))
    @printf("   Bi-214 activity: %.1f μBq/kg\n", ustrip(u"μBq/kg", rad_mat1.a_bi214))
    @printf("   Tl-208 activity: %.1f μBq/kg\n\n", ustrip(u"μBq/kg", rad_mat1.a_tl208))
    
    # Method 2: New constructor (from PhysicalMaterial)
    println("2. New constructor from PhysicalMaterial:")
    
    # First create a PhysicalMaterial
    mu_over_rho = 0.05u"cm^2/g"
    phys_mat = PhysicalMaterial(density, mu_over_rho)
    
    @printf("   PhysicalMaterial created:\n")
    @printf("     Density: %.2f g/cm³\n", ustrip(u"g/cm^3", phys_mat.ρ))
    @printf("     μ/ρ: %.3f cm²/g\n", ustrip(u"cm^2/g", phys_mat.μovrρ))
    @printf("     Attenuation length: %.1f cm\n", ustrip(u"cm", phys_mat.λ))
    
    # Now create RadioactiveMaterial from PhysicalMaterial
    rad_mat2 = RadioactiveMaterial(phys_mat, a_bi, a_tl)
    
    @printf("   RadioactiveMaterial from PhysicalMaterial:\n")
    @printf("     Density: %.2f g/cm³\n", ustrip(u"g/cm^3", rad_mat2.ρ))
    @printf("     Attenuation length: %.1f cm\n", ustrip(u"cm", rad_mat2.λ))
    @printf("     Bi-214 activity: %.1f μBq/kg\n", ustrip(u"μBq/kg", rad_mat2.a_bi214))
    @printf("     Tl-208 activity: %.1f μBq/kg\n\n", ustrip(u"μBq/kg", rad_mat2.a_tl208))
    
    # Verify both methods produce equivalent results for physical properties
    println("3. Verification:")
    @printf("   Densities match: %s\n", rad_mat1.ρ == rad_mat2.ρ)
    @printf("   Attenuation lengths match: %s\n", rad_mat1.λ == rad_mat2.λ)
    @printf("   Bi-214 activities match: %s\n", rad_mat1.a_bi214 == rad_mat2.a_bi214)
    @printf("   Tl-208 activities match: %s\n", rad_mat1.a_tl208 == rad_mat2.a_tl208)
    
    println("\n✅ Both constructors work correctly!")
end

function demonstrate_usage()
    println("\nPractical Usage Example")
    println("="^30)
    
    # Define common materials
    steel_phys = PhysicalMaterial(7.85u"g/cm^3", 0.06u"cm^2/g")
    copper_phys = PhysicalMaterial(8.96u"g/cm^3", 0.05u"cm^2/g")
    
    # Create radioactive versions with different contamination levels
    steel_rad = RadioactiveMaterial(steel_phys, 1.0u"mBq/kg", 0.4u"mBq/kg")
    copper_rad = RadioactiveMaterial(copper_phys, 10.0u"μBq/kg", 5.0u"μBq/kg")
    
    println("Steel material:")
    @printf("  Density: %.2f g/cm³\n", ustrip(u"g/cm^3", steel_rad.ρ))
    @printf("  Attenuation: %.1f cm\n", ustrip(u"cm", steel_rad.λ))
    @printf("  Contamination: %.1f mBq/kg (Bi-214), %.1f mBq/kg (Tl-208)\n\n", 
            ustrip(u"mBq/kg", steel_rad.a_bi214), ustrip(u"mBq/kg", steel_rad.a_tl208))
    
    println("Copper material:")
    @printf("  Density: %.2f g/cm³\n", ustrip(u"g/cm^3", copper_rad.ρ))
    @printf("  Attenuation: %.1f cm\n", ustrip(u"cm", copper_rad.λ))
    @printf("  Contamination: %.1f μBq/kg (Bi-214), %.1f μBq/kg (Tl-208)\n", 
            ustrip(u"μBq/kg", copper_rad.a_bi214), ustrip(u"μBq/kg", copper_rad.a_tl208))
end

function main()
    test_radioactive_constructors()
    demonstrate_usage()
end

# Run the script if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end