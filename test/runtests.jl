using SolarChem, Test
using StableRNGs

@testset "Load Data" begin include("load-data.jl") end
@testset "Data Management" begin include("data-mgmt.jl") end
@testset "Resamplers" begin include("resamplers.jl") end

