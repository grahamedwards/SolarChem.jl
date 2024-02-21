module SolarChem

import Random, DelimitedFiles
using VectorizedStatistics: vmean, vmedian
using NaNStatistics: nanmean, nanmedian, countnotnans

export periodictable
include("helpful-tuples.jl")

export loadastromatdata
include("load-data.jl")

export bootstrapmean, bootstrapresample
include("resamplers.jl")

export calcweights, estimateuncertainty, trimnans, pulltopic, pulltype, pullgroup, exclude
include("data-mgmt.jl")

end
