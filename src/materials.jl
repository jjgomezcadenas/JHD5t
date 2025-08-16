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

# Include materials data
include("MaterialsData.jl")


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

### Some materials:

# Helper function to create RadioactiveMaterial from material data
function create_material(name::String)
    if haskey(RHO, name) && haskey(MU, name)
        ρ = RHO[name]
        μovrρ = MU[name]
        μ = μovrρ * ρ
        λ = 1.0/μ
        
        # Get radioactivity data, default to zero if not found
        a_bi214 = get(BI214, name, 0.0*Bq/kg)
        a_tl208 = get(TL208, name, 0.0*Bq/kg)
        
        return RadioactiveMaterial(ρ, λ, a_bi214, a_tl208)
    else
        error("Material $name not found in data")
    end
end

# Pre-defined materials using data from MaterialsData.jl
fe316ti = create_material("Fe316Ti")  # Stainless steel 316Ti
copper = create_material("Cu")
ptfe = create_material("PTFE")
titanium = create_material("Ti")
lead = create_material("Pb")
    