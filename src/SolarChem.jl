module SolarChem

import  Random, Statistics, DelimitedFiles, Downloads, Serialization, Dates, # standard library
        Requires, # for data viz
        JSON3, # for handling data from astromat api
        Term, # to track astromat download progress
        Polyester, # for speed
        BonitoBook # for inspecting data

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

export cleanmetbase, cleanmetbase!, removemetbase, removemetbase!
include("metbase.jl")

include("mix.jl")
include("classic-stats.jl")
include("display.jl")

export inspect
"""
    inspect(filepath)

Generate and open a BonitoBook notebook environment to visually inspect Astromat data. Optionally provide a filepath to save the notebook to.

"""
function inspect(filepath::AbstractString="") 
    notebookname= "solarchem-visual-inspection"
    cpath = joinpath(ifelse(isempty(filepath),homedir(),filepath), string(notebookname,".md"))
    rm(cpath, force=true)

    rvi = read(joinpath(@__DIR__,"visual-inspection.md"), String)
    open(cpath, "w") do io
        write(io, rvi)
    end

    printstyled("\nBonitoBook notebook created at: ", color=:magenta)
    println("$cpath\n")
    printstyled("Warning: ",color=:yellow)
    print("this notebook will be overwritten by future executions of ")
    printstyled("inspect()", color=:cyan)
    println(". Save a copy to preserve any changes.\n\n")
    #printstyled("Notebook running at\n\n", bold=true)
    #println(string(x.protocol, x.url, ":", x.port, "/", notebookname))
    
    x = BonitoBook.book(cpath, openbrowser=false)
    
    printstyled("\n\nCopy-paste the link above into your browser to start inspecting!\n\n", color=:magenta)

end


function __init__()
    Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("visualization.jl")
    Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("visualization.jl") 
    Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" include("visualization.jl") 
end

end