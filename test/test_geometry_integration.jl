# Test integration between new geometry system and existing functionality
using JHD5t
using Test
using Unitful

@testset "Geometry Integration Tests" begin
    
    @testset "Real HD5t Detector Creation" begin
        # Test creating a realistic HD5t detector with the new geometry system
        
        # Create components as in the actual script
        bst_shell = CylinderShell(215.0u"cm", 220.0u"cm", 450.0u"cm")
        bst_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 1.0u"mBq/kg")
        bst = PhysicalCylindricalShell(bst_shell, bst_mat)
        
        bsl_shell = CylinderShell(200.0u"cm", 215.0u"cm", 400.0u"cm")
        bsl_mat = RadioactiveMaterial(8.96u"g/cm^3", 1.5u"cm", 10.0u"μBq/kg", 10.0u"μBq/kg")
        bsl = PhysicalCylindricalShell(bsl_shell, bsl_mat)
        
        est_cyl = Cylinder(220.0u"cm", 5.0u"cm")
        est_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 1.0u"mBq/kg")
        est = PhysicalCylinder(est_cyl, est_mat)
        
        esl_cyl = Cylinder(215.0u"cm", 15.0u"cm")
        esl_mat = RadioactiveMaterial(8.96u"g/cm^3", 1.5u"cm", 10.0u"μBq/kg", 10.0u"μBq/kg")
        esl = PhysicalCylinder(esl_cyl, esl_mat)
        
        gas_cyl = Cylinder(200.0u"cm", 400.0u"cm")
        gxe = GXe("rho_1520")
        gas_mat = RadioactiveMaterial(gxe)
        gas = PhysicalCylinder(gas_cyl, gas_mat)
        
        # Create vessel using backward-compatible constructor
        vessel = NextVessel(bst, bsl, est, esl, gas)
        
        # Test that envelope was created properly
        @test vessel.envelope isa Envelope
        @test vessel.envelope.shape isa Cylinder
        
        # Test envelope size is reasonable for the components
        env_bounds = envelope_bounds(vessel.envelope)
        @test ustrip(u"cm", env_bounds.xmax) >= 220.0  # At least as large as largest component
        @test ustrip(u"cm", env_bounds.zmax) >= 235.0  # At least as large as barrel + endcaps
        
        # Test all mass calculations still work
        total_mass_old = mass(bst) + mass(bsl) + mass(est) + mass(esl) + mass(gas)
        total_mass_new = mass(vessel)
        @test total_mass_new ≈ total_mass_old
        
        # Test individual mass functions
        @test mass_bst(vessel) ≈ mass(bst)
        @test mass_bsl(vessel) ≈ mass(bsl)
        @test mass_est(vessel) ≈ mass(est)
        @test mass_esl(vessel) ≈ mass(esl)
        @test mass_gas(vessel) ≈ mass(gas)
        
        # Test activity functions
        @test a_bi214_bst(vessel) == a_bi214(bst)
        @test a_tl208_bst(vessel) == a_tl208(bst)
        @test att_length_bst(vessel) == att_length(bst)
    end
    
    @testset "Custom Positioned Detector" begin
        # Test creating a detector with custom envelope and positioning
        
        # Create a large envelope
        envelope = CylindricalEnvelope(400.0u"cm", 800.0u"cm")
        
        # Create simple components
        shell1 = CylinderShell(50.0u"cm", 55.0u"cm", 100.0u"cm")
        mat1 = RadioactiveMaterial(1.0u"g/cm^3", 1.0u"cm", 0.0u"mBq/kg", 0.0u"mBq/kg")
        bst = PhysicalCylindricalShell(shell1, mat1)
        bsl = bst  # Reuse for simplicity
        
        cyl1 = Cylinder(60.0u"cm", 10.0u"cm")
        mat2 = RadioactiveMaterial(2.0u"g/cm^3", 1.0u"cm", 0.0u"mBq/kg", 0.0u"mBq/kg")
        est = PhysicalCylinder(cyl1, mat2)
        esl = est  # Reuse
        gas = est  # Reuse
        
        # Create custom positions
        bst_pos = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")      # Center
        bsl_pos = Position(10.0u"cm", 0.0u"cm", 0.0u"cm")     # Offset in X
        est_pos = Position(0.0u"cm", 0.0u"cm", 60.0u"cm")     # Positive Z
        esl_pos = Position(0.0u"cm", 0.0u"cm", -60.0u"cm")    # Negative Z
        gas_pos = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")      # Center
        
        # Create vessel with explicit positioning
        vessel = NextVessel(envelope, bst, bst_pos, bsl, bsl_pos, est, est_pos, esl, esl_pos, gas, gas_pos)
        
        # Test envelope preservation
        @test vessel.envelope == envelope
        
        # Test position preservation
        @test vessel.bst.position == bst_pos
        @test vessel.bsl.position == bsl_pos
        @test vessel.est.position == est_pos
        @test vessel.esl.position == esl_pos
        @test vessel.gas.position == gas_pos
        
        # Test volume names
        @test vessel.bst.name == "Barrel Structure"
        @test vessel.bsl.name == "Barrel Shield"
        @test vessel.est.name == "Endcap Structure"
        @test vessel.esl.name == "Endcap Shield"
        @test vessel.gas.name == "Gas Volume"
        
        # Test that physics calculations still work
        @test mass(vessel) > 0u"kg"
        @test mass_bst(vessel) == mass(bst)
    end
    
    @testset "Envelope Containment" begin
        # Test that all components fit within their envelope
        
        # Create a simple vessel
        shell = CylinderShell(10.0u"cm", 15.0u"cm", 20.0u"cm")
        mat = RadioactiveMaterial(1.0u"g/cm^3", 1.0u"cm", 0.0u"mBq/kg", 0.0u"mBq/kg")
        bst = PhysicalCylindricalShell(shell, mat)
        
        cyl = Cylinder(15.0u"cm", 5.0u"cm")
        est = PhysicalCylinder(cyl, mat)
        
        vessel = NextVessel(bst, bst, est, est, est)
        
        # Test that envelope bounds encompass all components
        env_bounds = envelope_bounds(vessel.envelope)
        
        # Test envelope is large enough for the largest component
        @test ustrip(u"cm", env_bounds.xmax) >= 15.0  # At least cylinder radius
        @test ustrip(u"cm", env_bounds.zmax) >= 15.0  # At least half the cylinder length + endcap
        
        # Test points within envelope
        center_point = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")
        @test is_inside(center_point, vessel.envelope)
        
        # Test points at component boundaries are within envelope
        edge_point = Position(14.0u"cm", 0.0u"cm", 0.0u"cm")
        @test is_inside(edge_point, vessel.envelope)
    end
    
    @testset "Volume Access by Name" begin
        # Test the volume lookup functionality
        
        shell = CylinderShell(5.0u"cm", 10.0u"cm", 15.0u"cm")
        mat = RadioactiveMaterial(1.0u"g/cm^3", 1.0u"cm", 0.0u"mBq/kg", 0.0u"mBq/kg")
        bst = PhysicalCylindricalShell(shell, mat)
        
        cyl = Cylinder(10.0u"cm", 5.0u"cm")
        est = PhysicalCylinder(cyl, mat)
        
        vessel = NextVessel(bst, bst, est, est, est)
        
        # Test getting volumes by name
        bst_vol = get_volume(vessel, "Barrel Structure")
        @test bst_vol == vessel.bst
        @test bst_vol.volume == bst
        
        bsl_vol = get_volume(vessel, "Barrel Shield")
        @test bsl_vol == vessel.bsl
        
        est_vol = get_volume(vessel, "Endcap Structure")
        @test est_vol == vessel.est
        
        esl_vol = get_volume(vessel, "Endcap Shield")
        @test esl_vol == vessel.esl
        
        gas_vol = get_volume(vessel, "Gas Volume")
        @test gas_vol == vessel.gas
        
        # Test getting non-existent volume
        missing_vol = get_volume(vessel, "Nonexistent Volume")
        @test missing_vol === nothing
        
        # Test placed_volumes function
        all_vols = placed_volumes(vessel)
        @test length(all_vols) == 5
        @test all_vols[1] == vessel.bst
        @test all_vols[2] == vessel.bsl
        @test all_vols[3] == vessel.est
        @test all_vols[4] == vessel.esl
        @test all_vols[5] == vessel.gas
    end
    
    @testset "Mass and Activity Consistency" begin
        # Test that mass and activity calculations are consistent between old and new systems
        
        # Create components with different materials
        bst_shell = CylinderShell(10.0u"cm", 15.0u"cm", 20.0u"cm")
        bst_mat = RadioactiveMaterial(7.0u"g/cm^3", 1.0u"cm", 5.0u"mBq/kg", 3.0u"mBq/kg")
        bst = PhysicalCylindricalShell(bst_shell, bst_mat)
        
        bsl_shell = CylinderShell(8.0u"cm", 10.0u"cm", 18.0u"cm")
        bsl_mat = RadioactiveMaterial(8.5u"g/cm^3", 2.0u"cm", 8.0u"μBq/kg", 12.0u"μBq/kg")
        bsl = PhysicalCylindricalShell(bsl_shell, bsl_mat)
        
        est_cyl = Cylinder(15.0u"cm", 3.0u"cm")
        est_mat = RadioactiveMaterial(7.2u"g/cm^3", 1.5u"cm", 2.0u"mBq/kg", 1.5u"mBq/kg")
        est = PhysicalCylinder(est_cyl, est_mat)
        
        esl_cyl = Cylinder(12.0u"cm", 8.0u"cm")
        esl_mat = RadioactiveMaterial(8.8u"g/cm^3", 2.5u"cm", 15.0u"μBq/kg", 20.0u"μBq/kg")
        esl = PhysicalCylinder(esl_cyl, esl_mat)
        
        gas_cyl = Cylinder(8.0u"cm", 18.0u"cm")
        gxe = GXe("rho_1520")
        gas_mat = RadioactiveMaterial(gxe)
        gas = PhysicalCylinder(gas_cyl, gas_mat)
        
        # Create vessel
        vessel = NextVessel(bst, bsl, est, esl, gas)
        
        # Test individual masses
        @test mass_bst(vessel) ≈ mass(bst)
        @test mass_bsl(vessel) ≈ mass(bsl)
        @test mass_est(vessel) ≈ mass(est)
        @test mass_esl(vessel) ≈ mass(esl)
        @test mass_gas(vessel) ≈ mass(gas)
        
        # Test total mass
        expected_total = mass(bst) + mass(bsl) + mass(est) + mass(esl) + mass(gas)
        @test mass(vessel) ≈ expected_total
        
        # Test activities
        @test a_bi214_bst(vessel) == a_bi214(bst)
        @test a_tl208_bst(vessel) == a_tl208(bst)
        @test a_bi214_bsl(vessel) == a_bi214(bsl)
        @test a_tl208_bsl(vessel) == a_tl208(bsl)
        @test a_bi214_est(vessel) == a_bi214(est)
        @test a_tl208_est(vessel) == a_tl208(est)
        @test a_bi214_esl(vessel) == a_bi214(esl)
        @test a_tl208_esl(vessel) == a_tl208(esl)
        
        # Test attenuation lengths
        @test att_length_bst(vessel) == att_length(bst)
        @test att_length_bsl(vessel) == att_length(bsl)
        @test att_length_est(vessel) == att_length(est)
        @test att_length_esl(vessel) == att_length(esl)
    end
    
    @testset "Envelope Volume Calculation" begin
        # Test that envelope volume calculations work correctly
        
        # Create envelope
        envelope = CylindricalEnvelope(50.0u"cm", 100.0u"cm")
        
        # Test envelope volume
        expected_volume = π * (50.0u"cm")^2 * 100.0u"cm"
        @test volume(envelope) ≈ expected_volume
        
        # Create a simple vessel and test envelope_volume function
        shell = CylinderShell(5.0u"cm", 10.0u"cm", 15.0u"cm")
        mat = RadioactiveMaterial(1.0u"g/cm^3", 1.0u"cm", 0.0u"mBq/kg", 0.0u"mBq/kg")
        bst = PhysicalCylindricalShell(shell, mat)
        
        cyl = Cylinder(10.0u"cm", 5.0u"cm")
        est = PhysicalCylinder(cyl, mat)
        
        vessel = NextVessel(bst, bst, est, est, est)
        
        # Test envelope_volume function
        @test envelope_volume(vessel) == volume(vessel.envelope)
        @test envelope_volume(vessel) > 0u"cm^3"
    end
end