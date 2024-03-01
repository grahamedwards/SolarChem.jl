using SolarChem, Test
using StableRNGs, Suppressor

@testset "Load Data" begin include("load-data.jl") end
@testset "Data Mgmt" begin include("data-mgmt.jl") end
@testset "Resamplers" begin include("resamplers.jl") end

