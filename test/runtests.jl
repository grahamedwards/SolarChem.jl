using SolarChem, Test
using StableRNGs

include("silence.jl")

mean(x) = SolarChem.StatsBase.mean(x)

@testset "Load Data" begin include("load-data.jl") end
@testset "Data Mgmt" begin include("data-mgmt.jl") end
@testset "Resamplers" begin include("resamplers.jl") end