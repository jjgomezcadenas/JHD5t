using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
using Interpolations

import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
	A, N, mol, mmol, V, L, M,
    Bq, mBq, μBq 

import PhysicalConstants.CODATA2018: N_A


A_BI214_316Ti   =    1.0    * mBq/kg
A_TL208_316Ti   =    0.4   * mBq/kg
A_BI214_CU_LIM  =   12      * μBq/kg
A_TL208_CU_LIM  =    1.4    * μBq/kg
A_BI214_CU_BEST =    3      * μBq/kg
A_TL208_CU_BEST =    1.4    * μBq/kg
A_BI214_PB      =  370      * μBq/kg
A_TL208_PB      =   73      * μBq/kg
A_BI214_Poly    =   62      * μBq/kg
A_TL208_Poly    =    8      * μBq/kg

rhoxe = Dict("rho_2020" => 124.3 * kg/m^3, 
               "rho_3020" => 203.35 * kg/m^3,
               "rho_1520" =>  89.9 * kg/m^3,
               "rho_1020" =>   58.0 * kg/m^3,
               "rho_0520" =>   30.0 * kg/m^3,
               "rho_0720" =>   40.0 * kg/m^3,
               "rho_0920" =>   50.0 * kg/m^3)

function rhoxegr(rhoxe_dict)
    return Dict(k => uconvert(g/cm^3, v) for (k, v) in rhoxe_dict)
end

rhoxe_gcm3 = rhoxegr(rhoxe)
"""
Simple representation of a physical material

# Fields
- `ρ::typeof(1.0g/cm^3)`          : density (ρ)
- `μovrρ::typeof(1.0cm^2/g)`      : μ/ρ
- `μ::Unitful.Length^-1`          : attenuation coefficient (μ) 
- `λ::Unitful.Length`             : attenuation length 

"""
struct PhysicalMaterial
    ρ::typeof(1.0g/cm^3) 
    μovrρ::typeof(1.0cm^2/g) 
	μ::typeof(1.0/cm)
	λ::Unitful.Length

    function PhysicalMaterial(ρ::typeof(1.0g/cm^3),  μovrρ::typeof(1.0cm^2/g))
        μ = μovrρ*ρ
        λ = 1.0/μ
		new(ρ, μovrρ, μ, λ)
	end
end

"""
Simple representation of a radioactive material

# Fields
- `m::PhysicalMaterial`          : a radioactive material is a physical material
- `a_bi214::typeof(1.0*Bq/kg)`   : activity of Bi-214
- `a_tl208::typeof(1.0*Bq/kg)`   : activity of Tl-208

"""
struct RadioactiveMaterial
    ρ::typeof(1.0g/cm^3)
    λ::Unitful.Length
    a_bi214::typeof(1.0*Bq/kg)
    a_tl208::typeof(1.0*Bq/kg)
    
    # Constructor from individual parameters
    RadioactiveMaterial(ρ, λ, a_bi214, a_tl208) = new(ρ, λ, a_bi214, a_tl208)
    
    # Constructor from PhysicalMaterial and radioactivity parameters
    function RadioactiveMaterial(phys_mat::PhysicalMaterial, a_bi214, a_tl208)
        new(phys_mat.ρ, phys_mat.λ, a_bi214, a_tl208)
    end
    
end

"""
Simple representation of xenon gas at P and T

# Fields
- `xe::PhysicalMaterial`         : xenon under conditions defined by dictionary

"""
struct GXe
    xe::PhysicalMaterial
    
    function GXe(ρ::String)
        
        if haskey(rhoxe_gcm3, ρ)
            xe = PhysicalMaterial(rhoxe_gcm3[ρ], 0.039 * cm^2/g)
        else 
            xe = PhysicalMaterial(rhoxe_gcm3["rho_1020"], 0.039 * cm^2/g)
        end
        new(xe)
    end 
end

# Constructor for RadioactiveMaterial from GXe with default zero radioactivity
function RadioactiveMaterial(gxe::GXe, a_bi214=0.0*Bq/kg, a_tl208=0.0*Bq/kg)
    RadioactiveMaterial(gxe.xe.ρ, gxe.xe.λ, a_bi214, a_tl208)
end