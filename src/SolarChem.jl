module SolarChem

import Random, DelimitedFiles, Requires
import Statistics
using Polyester: @batch # for speed

export periodictable, majors, heavysolar, lightsolar, mediumsolar, allsolar
include("helpful-tuples.jl")

export loadastromatdata, solartwins
include("load-data.jl")

export solartwins, sun
include("stars.jl")

export LoddersFegley1998
include("lodders.jl")

export bsmean, bsmean!, bsresample, bsresample!, bootstrapelements, bootstrapratios
include("resamplers.jl")

export calcweights, estimateuncertainty, trimnans, trimextremes, pulltopic, pulltype, pullgroup, exclude, excludeheated, countratios, countmeasurements, fraction2ratio, ratio2fraction, Composition, Fractions
include("data-mgmt.jl")

export ratiosummary
include("results.jl")

export solarlogmix, solarmixmetropolis
include("mix.jl")

include("display.jl")

function __init__()
    Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("visualization.jl")
    Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("visualization.jl") 
    Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" include("visualization.jl") 
end

end