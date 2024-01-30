module SolarChem

using DelimitedFiles: readdlm
using VectorizedStatistics: vmean, vmedian
using NaNStatistics: nanmean, nanmedian

# Write your package code here.
# include("load-data.jl")
# include(resamplers.jl"")

export calcweights, estimateuncertainty
include("data-mgmt.jl")

end
