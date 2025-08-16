"""
A Physical Cylinder shell includes:
1. A Cylinder shell shape
2. A radioactive material filling the shell
"""
struct PhysicalCylindricalShell
    shell::CylinderShell
    material::RadioactiveMaterial
end

inner_volume(c::PhysicalCylindricalShell) = inner_volume(c.shell)
shell_volume(c::PhysicalCylindricalShell) = shell_volume(c.shell)
inner_surface(c::PhysicalCylindricalShell) = inner_surface(c.shell)
outer_surface(c::PhysicalCylindricalShell) = outer_surface(c.shell)  
mass(c::PhysicalCylindricalShell) = shell_volume(c) * c.material.ρ
att_length(c::PhysicalCylindricalShell) = c.material.λ
a_bi214(c::PhysicalCylindricalShell) = c.material.a_bi214
a_tl208(c::PhysicalCylindricalShell) = c.material.a_tl208

"""
A Physical Cylinder includes:
1. A Cylinder shape
2. A radioactive material filling the cylinder
"""
struct PhysicalCylinder
    shell::Cylinder
    material::RadioactiveMaterial
end

volume(c::PhysicalCylinder) = volume(c.shell)
geometric_surface(c::PhysicalCylinder) = geometric_surface(c.shell)
endcap_surface(c::PhysicalCylinder) = endcap_surface(c.shell)
mass(c::PhysicalCylinder) = volume(c.shell) * c.material.ρ
att_length(c::PhysicalCylinder) = c.material.λ
a_bi214(c::PhysicalCylinder) = c.material.a_bi214
a_tl208(c::PhysicalCylinder) = c.material.a_tl208

"""
A NextVessel is made of:
1. An envelope (mother volume) that contains all components
2. A cylindrical shell shape defining the structure (e.g, steel)
3. A cylindrical shell shape defining the barrel shield (e.g, copper)
4. structure end-caps
5. shield end-caps
6. xenon
"""
struct NextVessel
    envelope::Envelope
    bst::PlacedVolume{PhysicalCylindricalShell}
    bsl::PlacedVolume{PhysicalCylindricalShell}
    est::PlacedVolume{PhysicalCylinder}
    esl::PlacedVolume{PhysicalCylinder}
    gas::PlacedVolume{PhysicalCylinder}
end

# Convenience constructor for backward compatibility
function NextVessel(bst::PhysicalCylindricalShell, bsl::PhysicalCylindricalShell,
                   est::PhysicalCylinder, esl::PhysicalCylinder, gas::PhysicalCylinder)
    
    # Create envelope that encompasses all components
    # Calculate required envelope size
    max_radius = max(bst.shell.Rout, est.shell.R)
    max_length = max(bst.shell.L, gas.shell.L) + 2 * est.shell.L  # Add endcap thickness
    
    envelope = CylindricalEnvelope(max_radius + 10.0u"cm", max_length + 20.0u"cm")
    
    # Create placed volumes at default positions (centered)
    bst_placed = PlacedVolume(bst, "Barrel Structure")
    bsl_placed = PlacedVolume(bsl, "Barrel Shield") 
    est_placed = PlacedVolume(est, "Endcap Structure")
    esl_placed = PlacedVolume(esl, "Endcap Shield")
    gas_placed = PlacedVolume(gas, "Gas Volume")
    
    return NextVessel(envelope, bst_placed, bsl_placed, est_placed, esl_placed, gas_placed)
end

# New constructor with explicit envelope and positioning
function NextVessel(envelope::Envelope,
                   bst::PhysicalCylindricalShell, bst_pos::Position,
                   bsl::PhysicalCylindricalShell, bsl_pos::Position,
                   est::PhysicalCylinder, est_pos::Position,
                   esl::PhysicalCylinder, esl_pos::Position,
                   gas::PhysicalCylinder, gas_pos::Position)
    
    bst_placed = PlacedVolume(bst, bst_pos, "Barrel Structure")
    bsl_placed = PlacedVolume(bsl, bsl_pos, "Barrel Shield")
    est_placed = PlacedVolume(est, est_pos, "Endcap Structure") 
    esl_placed = PlacedVolume(esl, esl_pos, "Endcap Shield")
    gas_placed = PlacedVolume(gas, gas_pos, "Gas Volume")
    
    return NextVessel(envelope, bst_placed, bsl_placed, est_placed, esl_placed, gas_placed)
end

# Helper functions for PlacedVolume
mass(pv::PlacedVolume) = mass(pv.volume)
a_bi214(pv::PlacedVolume) = a_bi214(pv.volume)
a_tl208(pv::PlacedVolume) = a_tl208(pv.volume)
att_length(pv::PlacedVolume) = att_length(pv.volume)

# Mass of individual components
mass_bst(n::NextVessel) = mass(n.bst)
mass_bsl(n::NextVessel) = mass(n.bsl)
mass_est(n::NextVessel) = mass(n.est)
mass_esl(n::NextVessel) = mass(n.esl)
mass_gas(n::NextVessel) = mass(n.gas)

# Total mass
mass(n::NextVessel) = mass_bst(n) + mass_bsl(n) + mass_est(n) + mass_esl(n) + mass_gas(n)

# Activity and attenuation properties for barrel structure
a_bi214_bst(n::NextVessel) = a_bi214(n.bst)
a_tl208_bst(n::NextVessel) = a_tl208(n.bst)
att_length_bst(n::NextVessel) = att_length(n.bst)

# Activity and attenuation properties for barrel shield
a_bi214_bsl(n::NextVessel) = a_bi214(n.bsl)
a_tl208_bsl(n::NextVessel) = a_tl208(n.bsl)
att_length_bsl(n::NextVessel) = att_length(n.bsl)

# Activity and attenuation properties for endcap structure
a_bi214_est(n::NextVessel) = a_bi214(n.est)
a_tl208_est(n::NextVessel) = a_tl208(n.est)
att_length_est(n::NextVessel) = att_length(n.est)

# Activity and attenuation properties for endcap shield
a_bi214_esl(n::NextVessel) = a_bi214(n.esl)
a_tl208_esl(n::NextVessel) = a_tl208(n.esl)
att_length_esl(n::NextVessel) = att_length(n.esl)

# Envelope and geometry utilities for NextVessel
envelope(n::NextVessel) = n.envelope
envelope_volume(n::NextVessel) = volume(n.envelope)
envelope_bounds(n::NextVessel) = envelope_bounds(n.envelope)

# Get all placed volumes as a collection
function placed_volumes(n::NextVessel)
    return [n.bst, n.bsl, n.est, n.esl, n.gas]
end

# Get volume by name
function get_volume(n::NextVessel, name::String)
    volumes = placed_volumes(n)
    for vol in volumes
        if vol.name == name
            return vol
        end
    end
    return nothing
end

     