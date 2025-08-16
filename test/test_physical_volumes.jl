# Test physical volumes functionality
using JHD5t
using Test
using Unitful

@testset "PhysicalVolume Tests" begin
    @testset "PhysicalCylindricalShell" begin
        # Create a shell and material
        shell = CylinderShell(9.0u"cm", 10.0u"cm", 100.0u"cm")
        mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 0.4u"mBq/kg")
        phys_shell = PhysicalCylindricalShell(shell, mat)
        
        # Test delegated methods
        @test inner_volume(phys_shell) == inner_volume(shell)
        @test shell_volume(phys_shell) == shell_volume(shell)
        @test inner_surface(phys_shell) == inner_surface(shell)
        @test outer_surface(phys_shell) == outer_surface(shell)
        
        # Test mass calculation
        expected_mass = shell_volume(shell) * mat.ρ
        @test mass(phys_shell) ≈ expected_mass
        
        # Test material properties
        @test att_length(phys_shell) == mat.λ
        @test a_bi214(phys_shell) == mat.a_bi214
        @test a_tl208(phys_shell) == mat.a_tl208
    end
    
    @testset "PhysicalCylinder" begin
        # Create a cylinder and material
        cyl = Cylinder(10.0u"cm", 50.0u"cm")
        mat = RadioactiveMaterial(8.96u"g/cm^3", 1.5u"cm", 3.0u"μBq/kg", 1.4u"μBq/kg")
        phys_cyl = PhysicalCylinder(cyl, mat)
        
        # Test delegated methods
        @test volume(phys_cyl) == volume(cyl)
        @test geometric_surface(phys_cyl) == geometric_surface(cyl)
        @test endcap_surface(phys_cyl) == endcap_surface(cyl)
        
        # Test mass calculation
        expected_mass = volume(cyl) * mat.ρ
        @test mass(phys_cyl) ≈ expected_mass
        
        # Test material properties
        @test att_length(phys_cyl) == mat.λ
        @test a_bi214(phys_cyl) == mat.a_bi214
        @test a_tl208(phys_cyl) == mat.a_tl208
    end
    
    @testset "NextVessel" begin
        # Create components for NextVessel
        bst_shell = CylinderShell(100.0u"cm", 101.0u"cm", 200.0u"cm")
        bst_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 0.4u"mBq/kg")
        bst = PhysicalCylindricalShell(bst_shell, bst_mat)
        
        bsl_shell = CylinderShell(101.0u"cm", 103.0u"cm", 200.0u"cm")
        bsl_mat = RadioactiveMaterial(8.96u"g/cm^3", 1.5u"cm", 3.0u"μBq/kg", 1.4u"μBq/kg")
        bsl = PhysicalCylindricalShell(bsl_shell, bsl_mat)
        
        est_cyl = Cylinder(100.0u"cm", 1.0u"cm")
        est_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 0.4u"mBq/kg")
        est = PhysicalCylinder(est_cyl, est_mat)
        
        esl_cyl = Cylinder(100.0u"cm", 2.0u"cm")
        esl_mat = RadioactiveMaterial(8.96u"g/cm^3", 1.5u"cm", 3.0u"μBq/kg", 1.4u"μBq/kg")
        esl = PhysicalCylinder(esl_cyl, esl_mat)
        
        gas_cyl = Cylinder(100.0u"cm", 200.0u"cm")
        gas_mat = RadioactiveMaterial(58.0u"kg/m^3", 25.6u"cm", 0.0u"Bq/kg", 0.0u"Bq/kg")
        gas = PhysicalCylinder(gas_cyl, gas_mat)
        
        vessel = NextVessel(bst, bsl, est, esl, gas)
        
        # Test individual mass functions
        @test mass_bst(vessel) == mass(bst)
        @test mass_bsl(vessel) == mass(bsl)
        @test mass_est(vessel) == mass(est)
        @test mass_esl(vessel) == mass(esl)
        @test mass_gas(vessel) == mass(gas)
        
        # Test total mass
        @test mass(vessel) == mass(bst) + mass(bsl) + mass(est) + mass(esl) + mass(gas)
        
        # Test barrel structure properties
        @test a_bi214_bst(vessel) == a_bi214(bst)
        @test a_tl208_bst(vessel) == a_tl208(bst)
        @test att_length_bst(vessel) == att_length(bst)
        
        # Test barrel shield properties
        @test a_bi214_bsl(vessel) == a_bi214(bsl)
        @test a_tl208_bsl(vessel) == a_tl208(bsl)
        @test att_length_bsl(vessel) == att_length(bsl)
        
        # Test endcap structure properties
        @test a_bi214_est(vessel) == a_bi214(est)
        @test a_tl208_est(vessel) == a_tl208(est)
        @test att_length_est(vessel) == att_length(est)
        
        # Test endcap shield properties
        @test a_bi214_esl(vessel) == a_bi214(esl)
        @test a_tl208_esl(vessel) == a_tl208(esl)
        @test att_length_esl(vessel) == att_length(esl)
    end
end