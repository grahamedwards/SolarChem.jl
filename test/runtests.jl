using SolarChem, Test
using StableRNGs

include("silence.jl")

mean(x) = SolarChem.Statistics.mean(x) 
std(x) = SolarChem.Statistics.std(x)

@testset "Load Data" include("load-data.jl")
@testset "Data Mgmt" include("data-mgmt.jl")
@testset "Resamplers" include("resamplers.jl")
@testset "Results" include("results.jl")
@testset "Stellar Data" include("stars.jl")
@testset "Lodd+Fegl,98" include("lodders.jl")