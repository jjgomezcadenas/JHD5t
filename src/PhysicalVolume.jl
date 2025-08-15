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
surface(c::PhysicalCylinder) = surface(c.shell)
endcap_surface(c::PhysicalCylinder) = endcap_surface(c.shell)
mass(c::PhysicalCylinder) = volume(c.shell) * c.material.ρ
att_length(c::PhysicalCylinder) = c.material.λ
a_bi214(c::PhysicalCylinder) = c.material.a_bi214
a_tl208(c::PhysicalCylinder) = c.material.a_tl208

"""
A NextVessel is made of:
1. A cylindrical shell shape defining the structure (e.g, steel)
2. A cylindrical shell shape defining the barrel shield (e.g, copper)
3. structure end-caps
4. shield end-caps
5. xenon
"""
struct NextVessel
    bst::PhysicalCylindricalShell
    bsl::PhysicalCylindricalShell
    est::PhysicalCylinder
    esl::PhysicalCylinder
    gas::PhysicalCylinder
end

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

     