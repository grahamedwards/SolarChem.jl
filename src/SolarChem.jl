module SolarChem

import Random, DelimitedFiles
using VectorizedStatistics: vmean, vmedian
using NaNStatistics: nanmean, nanmedian, countnotnans

export periodictable
include("helpful-tuples.jl")

export loadastromatdata
include("load-data.jl")

export bsmean, bsmean!, bsresample, bsresample!
include("resamplers.jl")

export calcweights, estimateuncertainty, trimnans, trimextremes, pulltopic, pulltype, pullgroup, exclude
include("data-mgmt.jl")

end
