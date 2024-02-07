module SolarChem

import Random, DelimitedFiles
using VectorizedStatistics: vmean, vmedian
using NaNStatistics: nanmean, nanmedian


include("helpful-tuples.jl")

include("load-data.jl")

export bootstrapmean, bootstrapresample
include("resamplers.jl")

export calcweights, estimateuncertainty, trimnans, pulltopic
include("data-mgmt.jl")

end
