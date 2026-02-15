"""
HD5t detector analysis functions

This module provides functions for reading and analyzing HD5t detector data,
including radioactivity budgets and detector component properties.
"""

using DataFrames
using CSV

"""
    read_detector_summary(filepath::String) -> DataFrame

Read the HD5t detector summary CSV file containing component masses and radioactivity data.

# Arguments
- `filepath::String`: Path to the hd5t_detector_summary.csv file

# Returns
- `DataFrame`: A DataFrame with the following columns:
  - Component: Name of the detector component
  - Mass_kg: Mass in kilograms
  - Bi214_Activity_mBq_kg: Bi-214 specific activity in mBq/kg
  - Tl208_Activity_mBq_kg: Tl-208 specific activity in mBq/kg
  - Total_Bi214_mBq: Total Bi-214 activity in mBq
  - Total_Tl208_mBq: Total Tl-208 activity in mBq

# Example
```julia
df = read_detector_summary("pluto/hd5t_detector_summary.csv")
```

# Notes
- The TOTAL row is included in the returned DataFrame
- NaN values in specific activity columns for TOTAL row are preserved
"""
function read_detector_summary(filepath::String)
    df = CSV.read(filepath, DataFrame)
    return df
end


"""
    get_component_data(df::DataFrame, component_name::String) -> DataFrameRow

Get the data for a specific detector component.

# Arguments
- `df::DataFrame`: DataFrame from read_detector_summary
- `component_name::String`: Name of the component (e.g., "Barrel_Shield_Cu")

# Returns
- `DataFrameRow`: Row containing all data for the specified component, or nothing if not found

# Example
```julia
df = read_detector_summary("pluto/hd5t_detector_summary.csv")
barrel_data = get_component_data(df, "Barrel_Shield_Cu")
println("Barrel mass: ", barrel_data.Mass_kg, " kg")
```
"""
function get_component_data(df::DataFrame, component_name::String)
    mask = df.Component .== component_name
    if sum(mask) == 0
        @warn "Component '$component_name' not found in detector summary"
        return nothing
    end
    return df[mask, :][1, :]
end


"""
    total_activity_bi214(df::DataFrame; exclude_total=true) -> Float64

Calculate total Bi-214 activity across all components.

# Arguments
- `df::DataFrame`: DataFrame from read_detector_summary
- `exclude_total::Bool`: If true, exclude the TOTAL row from calculation (default: true)

# Returns
- `Float64`: Total Bi-214 activity in mBq

# Example
```julia
df = read_detector_summary("pluto/hd5t_detector_summary.csv")
total_bi = total_activity_bi214(df)
println("Total Bi-214: ", total_bi, " mBq")
```
"""
function total_activity_bi214(df::DataFrame; exclude_total=true)
    data = exclude_total ? df[df.Component .!= "TOTAL", :] : df
    return sum(data.Total_Bi214_mBq)
end


"""
    total_activity_tl208(df::DataFrame; exclude_total=true) -> Float64

Calculate total Tl-208 activity across all components.

# Arguments
- `df::DataFrame`: DataFrame from read_detector_summary
- `exclude_total::Bool`: If true, exclude the TOTAL row from calculation (default: true)

# Returns
- `Float64`: Total Tl-208 activity in mBq

# Example
```julia
df = read_detector_summary("pluto/hd5t_detector_summary.csv")
total_tl = total_activity_tl208(df)
println("Total Tl-208: ", total_tl, " mBq")
```
"""
function total_activity_tl208(df::DataFrame; exclude_total=true)
    data = exclude_total ? df[df.Component .!= "TOTAL", :] : df
    return sum(data.Total_Tl208_mBq)
end


"""
    total_mass(df::DataFrame; exclude_total=true) -> Float64

Calculate total detector mass.

# Arguments
- `df::DataFrame`: DataFrame from read_detector_summary
- `exclude_total::Bool`: If true, exclude the TOTAL row from calculation (default: true)

# Returns
- `Float64`: Total mass in kg

# Example
```julia
df = read_detector_summary("pluto/hd5t_detector_summary.csv")
total_m = total_mass(df)
println("Total mass: ", total_m, " kg")
```
"""
function total_mass(df::DataFrame; exclude_total=true)
    data = exclude_total ? df[df.Component .!= "TOTAL", :] : df
    return sum(data.Mass_kg)
end


"""
    filter_by_material(df::DataFrame, material_suffix::String) -> DataFrame

Filter detector components by material type.

# Arguments
- `df::DataFrame`: DataFrame from read_detector_summary
- `material_suffix::String`: Material suffix (e.g., "Cu", "PTFE", "Fe316Ti", "Xe")

# Returns
- `DataFrame`: Filtered DataFrame containing only components with the specified material

# Example
```julia
df = read_detector_summary("pluto/hd5t_detector_summary.csv")
copper_components = filter_by_material(df, "Cu")
```
"""
function filter_by_material(df::DataFrame, material_suffix::String)
    mask = endswith.(df.Component, material_suffix)
    return df[mask, :]
end


"""
    activity_summary(df::DataFrame) -> DataFrame

Create a summary DataFrame showing total activities by material type.

# Arguments
- `df::DataFrame`: DataFrame from read_detector_summary

# Returns
- `DataFrame`: Summary with columns:
  - Material: Material type
  - Total_Bi214_mBq: Total Bi-214 activity for this material
  - Total_Tl208_mBq: Total Tl-208 activity for this material
  - Total_Mass_kg: Total mass for this material

# Example
```julia
df = read_detector_summary("pluto/hd5t_detector_summary.csv")
summary = activity_summary(df)
```
"""
function activity_summary(df::DataFrame)
    # Extract material types from component names
    materials = ["Cu", "PTFE", "Fe316Ti", "Xe"]

    results = DataFrame(
        Material = String[],
        Total_Bi214_mBq = Float64[],
        Total_Tl208_mBq = Float64[],
        Total_Mass_kg = Float64[]
    )

    for material in materials
        mat_df = filter_by_material(df, material)
        if nrow(mat_df) > 0
            push!(results, (
                Material = material,
                Total_Bi214_mBq = sum(mat_df.Total_Bi214_mBq),
                Total_Tl208_mBq = sum(mat_df.Total_Tl208_mBq),
                Total_Mass_kg = sum(mat_df.Mass_kg)
            ))
        end
    end

    return results
end
