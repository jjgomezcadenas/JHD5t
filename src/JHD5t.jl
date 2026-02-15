module JHD5t

using Unitful
using UnitfulEquivalences
using PhysicalConstants.CODATA2018
using Interpolations

include("shapes.jl")
include("materials.jl")
include("PhysicalVolume.jl")
include("histos.jl")
include("hd5t_functions.jl")
include("detector.jl")

export GeometricShape, Box, Cylinder, CylinderShell
export Position, Envelope, PlacedVolume, CylindricalEnvelope, BoxEnvelope
export volume, geometric_surface, endcap_surface, inner_volume, shell_volume
export inner_surface, inner_endcap_surface, outer_surface, outer_endcap_surface, thickness
export center_position, envelope_bounds, is_inside, ray_cylinder_intersection
export PhysicalMaterial, RadioactiveMaterial, GXe, create_material
export MU, RHO, BI214, TL208
export fe316ti, copper, ptfe, titanium, lead, hdpe, kapton
export PhysicalCylindricalShell, PhysicalCylinder, NextVessel
export mass, att_length, a_bi214, a_tl208
export mass_bst, mass_bsl, mass_est, mass_esl, mass_gas
export a_bi214_bst, a_tl208_bst, att_length_bst
export a_bi214_bsl, a_tl208_bsl, att_length_bsl
export a_bi214_est, a_tl208_est, att_length_est
export a_bi214_esl, a_tl208_esl, att_length_esl
export envelope, envelope_volume, placed_volumes, get_volume
export hist1d, hist2d, p1df, step_hist, in_range, get_histo1d, Histo1d
export scan_level, get_vals_from_sparse, centers, edges
export save_histo1d, load_histo1d
export save_histos, load_histos
export read_detector_summary, get_component_data
export total_activity_bi214, total_activity_tl208, total_mass
export filter_by_material, activity_summary
export hd5t_detector, analyze_hd5t_detector
export create_detector_summary_df, format_detector_summary_df, save_detector_summary
export plot_hd5t_projections, plot_hd5t_combined

end # module JHD5t
