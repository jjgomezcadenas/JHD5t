# Test package integration and exports
using JHD5t
using Test

@testset "Integration Tests" begin
    @testset "Basic Module Tests" begin
        @test typeof(JHD5t) == Module
    end
    
    @testset "Exported Functions" begin
        # Test that all exported functions are accessible
        
        # Geometry exports
        @test isdefined(JHD5t, :GeometricShape)
        @test isdefined(JHD5t, :Cylinder)
        @test isdefined(JHD5t, :CylinderShell)
        @test isdefined(JHD5t, :Box)
        
        # Material exports
        @test isdefined(JHD5t, :PhysicalMaterial)
        @test isdefined(JHD5t, :RadioactiveMaterial)
        @test isdefined(JHD5t, :GXe)
        
        # Physical volume exports
        @test isdefined(JHD5t, :PhysicalCylindricalShell)
        @test isdefined(JHD5t, :PhysicalCylinder)
        @test isdefined(JHD5t, :NextVessel)
        
        # Histogram exports
        @test isdefined(JHD5t, :hist1d)
        @test isdefined(JHD5t, :hist2d)
        @test isdefined(JHD5t, :Histo1d)
        @test isdefined(JHD5t, :save_histo1d)
        @test isdefined(JHD5t, :load_histo1d)
    end
end