module SolarChem

import  Random, Statistics, DelimitedFiles, Downloads, Serialization, Dates, # standard library
        Requires, # for data viz
        JSON3, # for handling data from astromat api
        Term, # to track astromat download progress
        Polyester # for speed

export periodictable, majors, heavysolar, lightsolar, mediumsolar, allsolar
include("helpful-tuples.jl")

export astromatdata, fastromat
include("get-data.jl")

export solartwins, sun
include("stars.jl")

export LoddersFegley1998
include("lodders.jl")

export bsmean, bsmean!, bsresample, bsresample!, bootstrapelements, bootstrapratios
include("resamplers.jl")

export calcweights, estimateuncertainty, trimnans, trimextremes, pulltopic, pulltype, pullgroup, exclude, excludeheated, countratios, countmeasurements, fraction2ratio, ratio2fraction, Composition, Fractions
include("data-mgmt.jl")

export ratiosummary, removefrom
include("results.jl")

export solarlogmix, solarmixmetropolis
include("mix.jl")

include("classic-stats.jl")

include("display.jl")

function __init__()
    Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("visualization.jl")
    Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("visualization.jl") 
    Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" include("visualization.jl") 

    Requires.@require Pluto="c3e4b0f8-55cb-11ea-2926-15256bba5781" begin
        """
            inspect()

        Visually inspect Astromat data in a Pluto notebook environment to identify specific data of interest.
        """
        inspect() = Pluto.run(notebook=joinpath(@__DIR__,"../visual-inspection.jl"))
        export inspect
    end
end

end