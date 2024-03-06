module SolarChem

import Random, DelimitedFiles, Requires
using VectorizedStatistics: vmean, vmedian
using NaNStatistics: nanmean, nanmedian, countnotnans

export periodictable
include("helpful-tuples.jl")

export loadastromatdata
include("load-data.jl")

export bsmean, bsmean!, bsresample, bsresample!, bootstrapelements, bootstrapratios
include("resamplers.jl")

export calcweights, estimateuncertainty, trimnans, trimextremes, pulltopic, pulltype, pullgroup, exclude, excludeheated, countratios, countmeasurements
include("data-mgmt.jl")

function __init__()
    Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("visualization.jl")
    Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("visualization.jl") 
end

end