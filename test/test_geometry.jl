# Test new geometry structures and functionality
using JHD5t
using Test
using Unitful

@testset "Geometry Tests" begin
    
    @testset "Position" begin
        # Test default constructor
        pos_default = Position()
        @test pos_default.x == 0.0u"cm"
        @test pos_default.y == 0.0u"cm"
        @test pos_default.z == 0.0u"cm"
        
        # Test constructor with coordinates
        pos = Position(10.0u"cm", 20.0u"cm", 30.0u"cm")
        @test pos.x == 10.0u"cm"
        @test pos.y == 20.0u"cm"
        @test pos.z == 30.0u"cm"
        
        # Test with different units
        pos_mm = Position(100.0u"mm", 200.0u"mm", 300.0u"mm")
        @test pos_mm.x == 100.0u"mm"
        @test pos_mm.y == 200.0u"mm"
        @test pos_mm.z == 300.0u"mm"
    end
    
    @testset "Box Volume" begin
        # Test Box volume calculation
        box = Box(0.0u"cm", 10.0u"cm", 0.0u"cm", 20.0u"cm", 0.0u"cm", 30.0u"cm")
        expected_volume = 10.0u"cm" * 20.0u"cm" * 30.0u"cm"
        @test volume(box) ≈ expected_volume
        
        # Test with negative coordinates
        box2 = Box(-5.0u"cm", 5.0u"cm", -10.0u"cm", 10.0u"cm", -15.0u"cm", 15.0u"cm")
        expected_volume2 = 10.0u"cm" * 20.0u"cm" * 30.0u"cm"
        @test volume(box2) ≈ expected_volume2
    end
    
    @testset "Envelope Creation" begin
        # Test CylindricalEnvelope with default center
        cyl_env = CylindricalEnvelope(50.0u"cm", 100.0u"cm")
        @test cyl_env.shape isa Cylinder
        @test cyl_env.shape.R == 50.0u"cm"
        @test cyl_env.shape.L == 100.0u"cm"
        @test cyl_env.center.x == 0.0u"cm"
        @test cyl_env.center.y == 0.0u"cm"
        @test cyl_env.center.z == 0.0u"cm"
        
        # Test CylindricalEnvelope with custom center
        center_pos = Position(10.0u"cm", 20.0u"cm", 30.0u"cm")
        cyl_env_offset = CylindricalEnvelope(25.0u"cm", 50.0u"cm", center_pos)
        @test cyl_env_offset.center.x == 10.0u"cm"
        @test cyl_env_offset.center.y == 20.0u"cm"
        @test cyl_env_offset.center.z == 30.0u"cm"
        
        # Test BoxEnvelope with default center
        box_env = BoxEnvelope(-10.0u"cm", 10.0u"cm", -20.0u"cm", 20.0u"cm", -30.0u"cm", 30.0u"cm")
        @test box_env.shape isa Box
        @test box_env.shape.xmin == -10.0u"cm"
        @test box_env.shape.xmax == 10.0u"cm"
        @test box_env.center.x == 0.0u"cm"
        
        # Test BoxEnvelope with custom center
        box_env_offset = BoxEnvelope(-5.0u"cm", 5.0u"cm", -10.0u"cm", 10.0u"cm", -15.0u"cm", 15.0u"cm", center_pos)
        @test box_env_offset.center == center_pos
    end
    
    @testset "Envelope Properties" begin
        # Test volume calculation
        cyl_env = CylindricalEnvelope(20.0u"cm", 40.0u"cm")
        expected_vol = π * (20.0u"cm")^2 * 40.0u"cm"
        @test volume(cyl_env) ≈ expected_vol
        
        # Test center position
        center_pos = Position(5.0u"cm", 10.0u"cm", 15.0u"cm")
        env_with_center = CylindricalEnvelope(30.0u"cm", 60.0u"cm", center_pos)
        @test center_position(env_with_center) == center_pos
    end
    
    @testset "Envelope Bounds" begin
        # Test cylindrical envelope bounds (centered at origin)
        cyl_env = CylindricalEnvelope(25.0u"cm", 50.0u"cm")
        bounds = envelope_bounds(cyl_env)
        @test bounds.xmin == -25.0u"cm"
        @test bounds.xmax == 25.0u"cm"
        @test bounds.ymin == -25.0u"cm"
        @test bounds.ymax == 25.0u"cm"
        @test bounds.zmin == -25.0u"cm"
        @test bounds.zmax == 25.0u"cm"
        
        # Test cylindrical envelope bounds (offset center)
        center_pos = Position(10.0u"cm", 20.0u"cm", 30.0u"cm")
        cyl_env_offset = CylindricalEnvelope(15.0u"cm", 40.0u"cm", center_pos)
        bounds_offset = envelope_bounds(cyl_env_offset)
        @test bounds_offset.xmin == -5.0u"cm"   # 10 - 15
        @test bounds_offset.xmax == 25.0u"cm"   # 10 + 15
        @test bounds_offset.ymin == 5.0u"cm"    # 20 - 15
        @test bounds_offset.ymax == 35.0u"cm"   # 20 + 15
        @test bounds_offset.zmin == 10.0u"cm"   # 30 - 20
        @test bounds_offset.zmax == 50.0u"cm"   # 30 + 20
        
        # Test box envelope bounds
        box_env = BoxEnvelope(-10.0u"cm", 10.0u"cm", -20.0u"cm", 20.0u"cm", -30.0u"cm", 30.0u"cm")
        box_bounds = envelope_bounds(box_env)
        @test box_bounds.xmin == -10.0u"cm"
        @test box_bounds.xmax == 10.0u"cm"
        @test box_bounds.ymin == -20.0u"cm"
        @test box_bounds.ymax == 20.0u"cm"
        @test box_bounds.zmin == -30.0u"cm"
        @test box_bounds.zmax == 30.0u"cm"
    end
    
    @testset "Point Inside Envelope" begin
        # Test point inside cylindrical envelope
        cyl_env = CylindricalEnvelope(25.0u"cm", 50.0u"cm")
        
        # Point at center should be inside
        center_point = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")
        @test is_inside(center_point, cyl_env) == true
        
        # Point within radius and length should be inside
        inside_point = Position(10.0u"cm", 10.0u"cm", 20.0u"cm")
        @test is_inside(inside_point, cyl_env) == true
        
        # Point outside radius should be outside
        outside_radius = Position(30.0u"cm", 0.0u"cm", 0.0u"cm")
        @test is_inside(outside_radius, cyl_env) == false
        
        # Point outside length should be outside
        outside_length = Position(0.0u"cm", 0.0u"cm", 30.0u"cm")
        @test is_inside(outside_length, cyl_env) == false
        
        # Test with box envelope
        box_env = BoxEnvelope(-10.0u"cm", 10.0u"cm", -20.0u"cm", 20.0u"cm", -30.0u"cm", 30.0u"cm")
        
        # Point inside box
        inside_box = Position(5.0u"cm", 15.0u"cm", 25.0u"cm")
        @test is_inside(inside_box, box_env) == true
        
        # Point outside box
        outside_box = Position(15.0u"cm", 5.0u"cm", 5.0u"cm")
        @test is_inside(outside_box, box_env) == false
    end
    
    @testset "PlacedVolume" begin
        # Create a test physical cylinder
        cyl = Cylinder(10.0u"cm", 20.0u"cm")
        mat = RadioactiveMaterial(1.0u"g/cm^3", 1.0u"cm", 0.0u"mBq/kg", 0.0u"mBq/kg")
        phys_cyl = PhysicalCylinder(cyl, mat)
        
        # Test PlacedVolume with default position
        placed_vol_default = PlacedVolume(phys_cyl, "Test Cylinder")
        @test placed_vol_default.volume == phys_cyl
        @test placed_vol_default.name == "Test Cylinder"
        @test placed_vol_default.position.x == 0.0u"cm"
        @test placed_vol_default.position.y == 0.0u"cm"
        @test placed_vol_default.position.z == 0.0u"cm"
        
        # Test PlacedVolume with custom position
        pos = Position(5.0u"cm", 10.0u"cm", 15.0u"cm")
        placed_vol = PlacedVolume(phys_cyl, pos, "Positioned Cylinder")
        @test placed_vol.volume == phys_cyl
        @test placed_vol.name == "Positioned Cylinder"
        @test placed_vol.position == pos
        
        # Test helper functions
        @test mass(placed_vol) == mass(phys_cyl)
        @test a_bi214(placed_vol) == a_bi214(phys_cyl)
        @test a_tl208(placed_vol) == a_tl208(phys_cyl)
        @test att_length(placed_vol) == att_length(phys_cyl)
    end
    
    @testset "NextVessel with Envelope - Backward Compatibility" begin
        # Create components the old way
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
        
        # Test old constructor - should create envelope automatically
        vessel = NextVessel(bst, bsl, est, esl, gas)
        
        # Verify envelope was created
        @test vessel.envelope isa Envelope
        @test vessel.envelope.shape isa Cylinder
        
        # Verify components are PlacedVolume objects
        @test vessel.bst isa PlacedVolume
        @test vessel.bsl isa PlacedVolume
        @test vessel.est isa PlacedVolume
        @test vessel.esl isa PlacedVolume
        @test vessel.gas isa PlacedVolume
        
        # Verify component names
        @test vessel.bst.name == "Barrel Structure"
        @test vessel.bsl.name == "Barrel Shield"
        @test vessel.est.name == "Endcap Structure"
        @test vessel.esl.name == "Endcap Shield"
        @test vessel.gas.name == "Gas Volume"
        
        # Verify all components are at default position (center)
        @test vessel.bst.position == Position()
        @test vessel.bsl.position == Position()
        @test vessel.est.position == Position()
        @test vessel.esl.position == Position()
        @test vessel.gas.position == Position()
        
        # Test envelope size (should be large enough to contain all components)
        env_bounds = envelope_bounds(vessel.envelope)
        @test env_bounds.xmax >= 220.0u"cm"  # At least as large as largest component
        @test env_bounds.zmax >= 235.0u"cm"  # At least as large as gas + endcaps
        
        # Test that old functions still work
        @test mass_bst(vessel) == mass(bst)
        @test mass_bsl(vessel) == mass(bsl)
        @test mass_est(vessel) == mass(est)
        @test mass_esl(vessel) == mass(esl)
        @test mass_gas(vessel) == mass(gas)
        
        # Test total mass
        expected_total = mass(bst) + mass(bsl) + mass(est) + mass(esl) + mass(gas)
        @test mass(vessel) ≈ expected_total
    end
    
    @testset "NextVessel with Explicit Envelope" begin
        # Create a custom envelope
        envelope = CylindricalEnvelope(300.0u"cm", 600.0u"cm")
        
        # Create components with custom positions
        bst_shell = CylinderShell(50.0u"cm", 55.0u"cm", 100.0u"cm")
        bst_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 1.0u"mBq/kg")
        bst = PhysicalCylindricalShell(bst_shell, bst_mat)
        bst_pos = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")
        
        # Minimal components for testing
        bsl = bst  # Reuse for simplicity
        bsl_pos = Position(10.0u"cm", 0.0u"cm", 0.0u"cm")
        
        est_cyl = Cylinder(60.0u"cm", 5.0u"cm")
        est_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 1.0u"mBq/kg")
        est = PhysicalCylinder(est_cyl, est_mat)
        est_pos = Position(0.0u"cm", 0.0u"cm", 52.5u"cm")
        
        esl = est  # Reuse for simplicity
        esl_pos = Position(0.0u"cm", 0.0u"cm", -52.5u"cm")
        
        gas_cyl = Cylinder(50.0u"cm", 100.0u"cm")
        gas_mat = RadioactiveMaterial(1.0u"g/cm^3", 1.0u"cm", 0.0u"mBq/kg", 0.0u"mBq/kg")
        gas = PhysicalCylinder(gas_cyl, gas_mat)
        gas_pos = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")
        
        # Create vessel with explicit envelope and positions
        vessel = NextVessel(envelope, bst, bst_pos, bsl, bsl_pos, est, est_pos, esl, esl_pos, gas, gas_pos)
        
        # Test that envelope is preserved
        @test vessel.envelope == envelope
        
        # Test that positions are preserved
        @test vessel.bst.position == bst_pos
        @test vessel.bsl.position == bsl_pos
        @test vessel.est.position == est_pos
        @test vessel.esl.position == esl_pos
        @test vessel.gas.position == gas_pos
    end
    
    @testset "NextVessel Utility Functions" begin
        # Create test components with correct types
        bst_shell = CylinderShell(10.0u"cm", 15.0u"cm", 20.0u"cm")
        bst_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 1.0u"mBq/kg")
        bst = PhysicalCylindricalShell(bst_shell, bst_mat)
        bsl = bst  # Reuse shell for bsl
        
        # Create cylinder components for endcaps and gas
        cyl = Cylinder(15.0u"cm", 5.0u"cm")
        cyl_mat = RadioactiveMaterial(7.85u"g/cm^3", 2.0u"cm", 1.0u"mBq/kg", 1.0u"mBq/kg")
        est = PhysicalCylinder(cyl, cyl_mat)
        esl = est  # Reuse cylinder for esl
        gas = est  # Reuse cylinder for gas
        
        vessel = NextVessel(bst, bsl, est, esl, gas)
        
        # Test envelope accessors
        @test envelope(vessel) == vessel.envelope
        @test envelope_volume(vessel) == volume(vessel.envelope)
        @test envelope_bounds(vessel) == envelope_bounds(vessel.envelope)
        
        # Test placed_volumes function
        volumes = placed_volumes(vessel)
        @test length(volumes) == 5
        @test volumes[1] == vessel.bst
        @test volumes[2] == vessel.bsl
        @test volumes[3] == vessel.est
        @test volumes[4] == vessel.esl
        @test volumes[5] == vessel.gas
        
        # Test get_volume function
        @test get_volume(vessel, "Barrel Structure") == vessel.bst
        @test get_volume(vessel, "Barrel Shield") == vessel.bsl
        @test get_volume(vessel, "Endcap Structure") == vessel.est
        @test get_volume(vessel, "Endcap Shield") == vessel.esl
        @test get_volume(vessel, "Gas Volume") == vessel.gas
        @test get_volume(vessel, "Nonexistent") === nothing
    end
    
    @testset "Ray-Cylinder Intersection Placeholder" begin
        # Test that the placeholder function exists and returns empty array
        ray_origin = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")
        ray_direction = Position(1.0u"cm", 0.0u"cm", 0.0u"cm")
        cylinder = Cylinder(10.0u"cm", 20.0u"cm")
        cylinder_center = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")
        
        intersections = ray_cylinder_intersection(ray_origin, ray_direction, cylinder, cylinder_center)
        @test intersections isa Vector{Float64}
        @test length(intersections) == 0  # Placeholder returns empty array
    end
end