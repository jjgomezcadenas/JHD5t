# Test materials functionality
using JHD5t
using Test
using Unitful

@testset "Materials Tests" begin
    @testset "PhysicalMaterial" begin
        # Test material creation
        density = 7.85u"g/cm^3"  # Steel density
        mu_over_rho = 0.05u"cm^2/g"
        mat = PhysicalMaterial(density, mu_over_rho)
        
        @test mat.ρ == density
        @test mat.μovrρ == mu_over_rho
        @test mat.μ ≈ mu_over_rho * density
        @test mat.λ ≈ 1.0 / (mu_over_rho * density)
    end
    
    @testset "RadioactiveMaterial" begin
        density = 8.96u"g/cm^3"  # Copper density
        att_length = 1.0u"cm"
        a_bi = 3.0u"μBq/kg"
        a_tl = 1.4u"μBq/kg"
        
        # Test constructor from individual parameters
        rad_mat = RadioactiveMaterial(density, att_length, a_bi, a_tl)
        
        @test rad_mat.ρ == density
        @test rad_mat.λ == att_length
        @test rad_mat.a_bi214 == a_bi
        @test rad_mat.a_tl208 == a_tl
        
        # Test constructor from PhysicalMaterial
        phys_mat = PhysicalMaterial(density, 0.05u"cm^2/g")
        rad_mat2 = RadioactiveMaterial(phys_mat, a_bi, a_tl)
        
        @test rad_mat2.ρ == phys_mat.ρ
        @test rad_mat2.λ == phys_mat.λ
        @test rad_mat2.a_bi214 == a_bi
        @test rad_mat2.a_tl208 == a_tl
        
        # Test that both constructors give equivalent results when given same parameters
        rad_mat3 = RadioactiveMaterial(phys_mat.ρ, phys_mat.λ, a_bi, a_tl)
        @test rad_mat2.ρ == rad_mat3.ρ
        @test rad_mat2.λ == rad_mat3.λ
        @test rad_mat2.a_bi214 == rad_mat3.a_bi214
        @test rad_mat2.a_tl208 == rad_mat3.a_tl208
    end
    
    @testset "GXe" begin
        # Test default xenon density
        gxe = GXe("rho_1020")
        @test gxe.xe.ρ == 0.058u"g/cm^3"
        @test gxe.xe.μovrρ == 0.039u"cm^2/g"
        
        # Test invalid key defaults to rho_1020
        gxe_default = GXe("invalid_key")
        @test gxe_default.xe.ρ == 0.058u"g/cm^3"
        
        # Test different density options
        gxe_1520 = GXe("rho_1520")
        @test gxe_1520.xe.ρ ≈ 0.0899u"g/cm^3"  # 89.9 kg/m³ = 0.0899 g/cm³
        
        # Test that GXe can be converted to RadioactiveMaterial
        rad_mat_from_gxe = RadioactiveMaterial(gxe)
        @test rad_mat_from_gxe.ρ == gxe.xe.ρ
        @test rad_mat_from_gxe.λ == gxe.xe.λ
        @test rad_mat_from_gxe.a_bi214 == 0.0u"Bq/kg"
        @test rad_mat_from_gxe.a_tl208 == 0.0u"Bq/kg"
        
        # Test GXe to RadioactiveMaterial with custom activities
        custom_rad_mat = RadioactiveMaterial(gxe, 1.0u"μBq/kg", 0.5u"μBq/kg")
        @test custom_rad_mat.ρ == gxe.xe.ρ
        @test custom_rad_mat.λ == gxe.xe.λ
        @test custom_rad_mat.a_bi214 == 1.0u"μBq/kg"
        @test custom_rad_mat.a_tl208 == 0.5u"μBq/kg"
    end
end