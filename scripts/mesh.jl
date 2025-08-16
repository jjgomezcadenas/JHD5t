#!/usr/bin/env julia

# Hex-mesh utilities for a perforated stainless-steel disk
# Assumes: hexagonal holes with across-flats inner "diameter" D,
#          strut (wire) width w, sheet thickness t_mesh, radius R, density rho.

using Printf

"""
    solid_fraction_hex(D, w)

Solid metal area fraction of a hex mesh made by etching hexagonal holes.
`D` = across-flats inner hex size (m). `w` = strut width (m).

Model: unit cell across-flats pitch = D + w.
f = 1 - (D/(D+w))^2
"""
solid_fraction_hex(D, w) = 1 - (D / (D + w))^2

"""
    mesh_mass(R, t_mesh, rho, D, w)

Mass (kg) of a circular hex mesh of radius `R` (m), thickness `t_mesh` (m),
material density `rho` (kg/m^3), hex inner size `D` (m), strut width `w` (m).
"""
function mesh_mass(R, t_mesh, rho, D, w)
    f = solid_fraction_hex(D, w)
    A = π * R^2
    V = f * A * t_mesh
    return rho * V
end

"""
    equivalent_solid_thickness(t_mesh, D, w)

Thickness (m) of a *solid* disk (same radius & material) that has the same mass
as the hex mesh of thickness `t_mesh`. Depends only on the solid-area fraction.
"""
equivalent_solid_thickness(t_mesh, D, w) = solid_fraction_hex(D, w) * t_mesh

# -------------------------------
# Example with your numbers
# -------------------------------
R      = 2000.0 / 1000       # 2000 mm → 2.0 m
t_mesh = 127e-6              # 127 μm
D      = 5e-3                # 5 mm inner hex across-flats
w      = 0.127e-3            # 127 μm strut width
rho    = 8000.0              # kg/m^3 (stainless steel)

f       = solid_fraction_hex(D, w)
m_mesh  = mesh_mass(R, t_mesh, rho, D, w)
t_solid = equivalent_solid_thickness(t_mesh, D, w)

println("\n=== Hex Mesh Calculations ===")
println("Parameters:")
@printf("  Radius: %.1f mm (%.3f m)\n", R*1000, R)
@printf("  Mesh thickness: %.1f μm (%.3e m)\n", t_mesh*1e6, t_mesh)
@printf("  Hex hole size (across-flats): %.1f mm\n", D*1000)
@printf("  Strut width: %.1f μm\n", w*1e6)
@printf("  Material density: %.0f kg/m³\n", rho)

println("\nResults:")
@printf("  Solid fraction: %.3f%% metal\n", f*100)
@printf("  Mesh mass: %.3f kg\n", m_mesh)
@printf("  Equivalent solid thickness: %.2f μm\n", t_solid*1e6)