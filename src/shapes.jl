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

"""
    An abstract class representing a geometrical shape.
"""
abstract type GeometricShape end

"""
    Represents a box of dimensions (xmin,xmax), (ymin, ymax), (zmin, zmax)
"""
struct Box <: GeometricShape
    xmin::Unitful.Length 
    xmax::Unitful.Length 
    ymin::Unitful.Length 
    ymax::Unitful.Length 
    zmin::Unitful.Length 
    zmax::Unitful.Length 
end


"""
    Represents a cylinder of radius R and length L
"""
struct Cylinder <: GeometricShape
    R::Unitful.Length 
    L::Unitful.Length 
end



"""
    Represents a cylinder shell of internal radius Rin
    external radius Rout and length L
"""
struct CylinderShell <: GeometricShape
    Rin::Unitful.Length 
    Rout::Unitful.Length 
    L::Unitful.Length 
end


volume(c::Cylinder) = π * c.R^2 * c.L

geometric_surface(c::Cylinder) = 2π * c.R * c.L

endcap_surface(c::Cylinder) = π * c.R^2

inner_volume(c::CylinderShell) = π * c.Rin^2 * c.L

shell_volume(c::CylinderShell) = π * (c.Rout^2 - c.Rin^2) * c.L
    
inner_surface(c::CylinderShell) = 2π * c.Rin * c.L  

inner_endcap_surface(c::CylinderShell) = π * c.Rin^2

outer_surface(c::CylinderShell) = 2π * c.Rout * c.L

outer_endcap_surface(c::CylinderShell) = π * c.Rout^2

thickness(c::CylinderShell) = c.Rout - c.Rin

# Volume functions for Box
volume(b::Box) = (b.xmax - b.xmin) * (b.ymax - b.ymin) * (b.zmax - b.zmin)

"""
    Position represents a 3D position in space with optional rotation
"""
struct Position
    x::Unitful.Length
    y::Unitful.Length  
    z::Unitful.Length
    # Future: could add rotation angles (θx, θy, θz)
end

Position() = Position(0.0u"cm", 0.0u"cm", 0.0u"cm")

"""
    Envelope represents a mother volume that contains all detector components
    Currently supports cylindrical envelopes
"""
struct Envelope
    shape::GeometricShape
    center::Position
end

# Create cylindrical envelope
function CylindricalEnvelope(radius::Unitful.Length, length::Unitful.Length, 
                           center::Position = Position())
    return Envelope(Cylinder(radius, length), center)
end

# Create box envelope  
function BoxEnvelope(xmin::Unitful.Length, xmax::Unitful.Length,
                    ymin::Unitful.Length, ymax::Unitful.Length,
                    zmin::Unitful.Length, zmax::Unitful.Length,
                    center::Position = Position())
    return Envelope(Box(xmin, xmax, ymin, ymax, zmin, zmax), center)
end

# Envelope properties
volume(env::Envelope) = volume(env.shape)
center_position(env::Envelope) = env.center

"""
    PlacedVolume represents a physical volume at a specific position within an envelope
"""
struct PlacedVolume{T}
    volume::T  # Can be PhysicalCylinder, PhysicalCylindricalShell, etc.
    position::Position
    name::String
end

# Convenience constructor with default position
PlacedVolume(volume, name::String) = PlacedVolume(volume, Position(), name)

# Geometric queries
function is_inside(point::Position, envelope::Envelope)
    # Check if point is inside the envelope
    if envelope.shape isa Cylinder
        cyl = envelope.shape
        # Distance from center axis
        r = sqrt((point.x - envelope.center.x)^2 + (point.y - envelope.center.y)^2)
        # Z position relative to envelope center
        z_rel = abs(point.z - envelope.center.z)
        return r <= cyl.R && z_rel <= cyl.L/2
    elseif envelope.shape isa Box
        box = envelope.shape
        return (box.xmin <= point.x <= box.xmax &&
                box.ymin <= point.y <= box.ymax &&
                box.zmin <= point.z <= box.zmax)
    end
    return false
end

function envelope_bounds(envelope::Envelope)
    """Get the bounding box of an envelope"""
    if envelope.shape isa Cylinder
        cyl = envelope.shape
        cx, cy, cz = envelope.center.x, envelope.center.y, envelope.center.z
        return (
            xmin = cx - cyl.R, xmax = cx + cyl.R,
            ymin = cy - cyl.R, ymax = cy + cyl.R, 
            zmin = cz - cyl.L/2, zmax = cz + cyl.L/2
        )
    elseif envelope.shape isa Box
        box = envelope.shape
        cx, cy, cz = envelope.center.x, envelope.center.y, envelope.center.z
        return (
            xmin = cx + box.xmin, xmax = cx + box.xmax,
            ymin = cy + box.ymin, ymax = cy + box.ymax,
            zmin = cz + box.zmin, zmax = cz + box.zmax
        )
    end
end

# Ray-geometry intersection utilities (for future ray tracing)
function ray_cylinder_intersection(ray_origin::Position, ray_direction::Position, 
                                 cylinder::Cylinder, cylinder_center::Position)
    """
    Find intersection points of a ray with a cylinder
    Returns distances along ray to intersection points, or empty if no intersection
    """
    # This is a placeholder for future ray tracing implementation
    # Would implement the mathematical ray-cylinder intersection algorithm
    return Float64[]
end

