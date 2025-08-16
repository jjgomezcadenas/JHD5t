# Main test runner for JHD5t package
using JHD5t
using Test

@testset "JHD5t.jl" begin
    # Run all test suites
    include("test_integration.jl")
    include("test_shapes.jl")
    include("test_materials.jl")
    include("test_physical_volumes.jl")
    include("test_histograms.jl")
    include("test_geometry.jl")
    include("test_geometry_integration.jl")
end