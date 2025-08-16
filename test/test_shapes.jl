# Test geometric shapes functionality
using JHD5t
using Test
using Unitful

@testset "Shapes Tests" begin
    @testset "Cylinder" begin
        cyl = Cylinder(10.0u"cm", 20.0u"cm")
        @test cyl.R == 10.0u"cm"
        @test cyl.L == 20.0u"cm"
        
        # Test volume calculation
        expected_volume = π * (10.0u"cm")^2 * 20.0u"cm"
        @test volume(cyl) ≈ expected_volume
        
        # Test surface calculation (lateral surface)
        expected_surface = 2π * 10.0u"cm" * 20.0u"cm"
        @test geometric_surface(cyl) ≈ expected_surface
        
        # Test endcap surface
        expected_endcap = π * (10.0u"cm")^2
        @test endcap_surface(cyl) ≈ expected_endcap
    end
    
    @testset "CylinderShell" begin
        shell = CylinderShell(8.0u"cm", 10.0u"cm", 20.0u"cm")
        @test shell.Rin == 8.0u"cm"
        @test shell.Rout == 10.0u"cm"
        @test shell.L == 20.0u"cm"
        
        # Test inner volume
        expected_inner_vol = π * (8.0u"cm")^2 * 20.0u"cm"
        @test inner_volume(shell) ≈ expected_inner_vol
        
        # Test shell volume
        expected_shell_vol = π * ((10.0u"cm")^2 - (8.0u"cm")^2) * 20.0u"cm"
        @test shell_volume(shell) ≈ expected_shell_vol
        
        # Test inner surface
        expected_inner_surf = 2π * 8.0u"cm" * 20.0u"cm"
        @test inner_surface(shell) ≈ expected_inner_surf
        
        # Test outer surface
        expected_outer_surf = 2π * 10.0u"cm" * 20.0u"cm"
        @test outer_surface(shell) ≈ expected_outer_surf
        
        # Test thickness
        @test thickness(shell) == 2.0u"cm"
        
        # Test endcap surfaces
        @test inner_endcap_surface(shell) ≈ π * (8.0u"cm")^2
        @test outer_endcap_surface(shell) ≈ π * (10.0u"cm")^2
    end
    
    @testset "Box" begin
        box = Box(0.0u"cm", 10.0u"cm", 0.0u"cm", 20.0u"cm", 0.0u"cm", 30.0u"cm")
        @test box.xmin == 0.0u"cm"
        @test box.xmax == 10.0u"cm"
        @test box.ymin == 0.0u"cm"
        @test box.ymax == 20.0u"cm"
        @test box.zmin == 0.0u"cm"
        @test box.zmax == 30.0u"cm"
        
        # Test volume calculation
        expected_volume = 10.0u"cm" * 20.0u"cm" * 30.0u"cm"
        @test volume(box) ≈ expected_volume
    end
end