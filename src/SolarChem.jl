module SolarChem

import Random
using DelimitedFiles: readdlm
using VectorizedStatistics: vmean, vmedian
using NaNStatistics: nanmean, nanmedian

# Write your package code here.
include("load-data.jl")

export bootstrapmean, bootstrapresample
include("resamplers.jl")

export calcweights, estimateuncertainty, trimnans, pulltopic
include("data-mgmt.jl")

end