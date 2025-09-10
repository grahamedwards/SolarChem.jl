
"""

```julia
    astromatdata(; url, datasets_per_request=1000, save=true, showskipped=false)
```

Get elemental geochemical data (only) from the AstroMat database for all meteorites, returned as a NamedTuple. The default `url` string is:

"https://api.astromat.org/v4/search/results?analysisTypes=rock%3A%3A[WHOLE+ROCK]&taxons=METEORITE&variables=MAJ||REE||TE||VO"

This requests whole rock major (MJ), volatile (VO), trace (TE), and rare earth element (REE) data for all meteorites in the database. You may provide your own customized API string, but this is not recommended unless you are familiar with the AstroMat API.

By default, the function requests 1000 `datasets_per_request`, which is a roughly optimized amount.

When `save` is set to `true` or provided with a filepath, the data is serlialized to the provided filepath or to a file automatically labeled `astromat-YYYY-MM-DD.jls` (reflecting the current date). If `false` is provided, the output is not saved.

When `showskipped` is `true`, returns a tuple of the data from astromat, as well as a list of excluded variable names.

"""
function astromatdata(; url::String="https://api.astromat.org/v4/search/results?analysisTypes=rock::[WHOLE+ROCK]&taxons=METEORITE&variables=MAJ||REE||TE||VO", datasets_per_request::Int = 1000, save::Union{String,Bool}=true, showskipped::Bool=false)

    
    url = string(url,"&size=$datasets_per_request")

    save, filesavelocation = if save isa AbstractString 
        fp, _ = splitdir(save)
        @assert isdir(filepath) "Invalid filepath provided: '$fp/' does not exist"
        true, save
    elseif save
        true, joinpath(pwd(),"astromat-$(Dates.today()).jls")
    else
        false, nothing
    end

    numlookup = let
        tnn = DelimitedFiles.readdlm(joinpath(@__DIR__,"..","data","astromat","taxon-num-name.txt"),'\t', skipstart=1)
        nums = Vector{Union{String,Int}}(undef, size(tnn)[1])
        rawtaxon = Vector{String}(undef, size(tnn)[1])
        tps = Vector{Vector{Int}}(undef,size(tnn)[1])
        grps = Vector{Vector{Symbol}}(undef,size(tnn)[1])
        @inbounds for i = axes(tnn)[1]
            num, grp, tp, txn = @view tnn[i,:]
            nums[i] = Int(num) #isa Number ? Int(num) : String(num)
            grps[i] = Symbol.(split(grp, ','))
            tps[i] = digits(ifelse(isempty(tp),0, tp))
            rawtaxon[i] = string(txn)
        end 

        Dict( nums .=> [[grps[i], tps[i], rawtaxon[i]] for i = eachindex(nums)])
    end


    ellookup = Dict(string.(periodictable()) .=> periodictable())
    ellookup["FeTOT"] = :Fe

    oxidelookup = Dict("CaO" => :Ca, "Na2O" => :Na, "NiO" => :Ni, "CoO" => :Co, "MgO" => :Mg, "FeO" => :Fe,  "Fe2O3" => :Fe3, "SiO2" => :Si, "Al2O3" => :Al, "MnO" => :Mn, "Cr2O3" => :Cr, "TiO2" => :Ti, "K2O" => :K , "P2O5" => :P)


    # https://astromat.org/synthesis?analysisTypes=rock::[WHOLE+ROCK]&taxons=METEORITE&variables=MAJ||REE||TE||VO


    io = IOBuffer()
    jsn = JSON3.read( String( take!( Downloads.download(url, io) ) ) )# formerly String(HTTP.get(url).body)

    n = jsn[:count]

    name = Vector{String}(undef, n)
    comment, citation, dataset = similar(name), similar(name), similar(name)
    group = Vector{Vector{Symbol}}(undef,n)
    type = Vector{Vector{Int}}(undef,n)
    taxon = similar(name)

    meas = NamedTuple{periodictable()}(fill(NaN,n) for i = eachindex(periodictable()))
    uncs = deepcopy(meas)

    ignored = []

    ind = 0 # running index of output dataset

    pbar = Term.Progress.ProgressBar(;columns=:detailed)
    job = Term.Progress.addjob!(pbar; N=ceil(Int,n/datasets_per_request))
    Term.Progress.start!(pbar)
    Term.Progress.render(pbar)
    while !isnothing(jsn[:afterKey])
        @inbounds for i = eachindex(jsn[:data])
            ind += 1
            temptest = jsn[:data][i][:temperature]
            heated = jsn[:data][i][:heated]
            isnothing(temptest) || @warn "Non-nothing temperature of $temptest"
            isnothing(heated) || @warn "heated value of $heated"


            name[ind] = jsn[:data][i][:sample][:sampleCode] # sampleCode is the official, normalized name as it is recognized in MetBull that was linked to the sampleName (the name that was given in the paper) by the curator.
            comment[ind] = ifelse(isnothing(jsn[:data][i][:analysisComment]), "", jsn[:data][i][:analysisComment])
            citation[ind] = string("https://astromat.org/synthesis/citation/",jsn[:data][i][:citation][:id])
            dataset[ind] = jsn[:data][i][:dataset][:datasetTitle]
            group[ind], type[ind], taxon[ind] = numlookup[jsn[:data][i][:sample][:taxon][1][:num]]

            @inbounds for j = eachindex(jsn[:data][i][:variables])
                variable = jsn[:data][i][:variables][j][:variable]
                value = float(jsn[:data][i][:variables][j][:valueMeas])
                unit = jsn[:data][i][:variables][j][:unit]

                stdevUnit = jsn[:data][i][:variables][j][:stdevUnit]
                
                stdevType = jsn[:data][i][:variables][j][:stdevType]
                stdevType = ifelse(isnothing(stdevType), "", stdevType)

                stdevValue = jsn[:data][i][:variables][j][:stdevValue]
                try 
                    stdevValue =  parse(Float64, stdevValue)
                catch
                    stdevValue = NaN 
                end

                value *= SolarChem.unitconversionfactor(unit; io=devnull)

                stdevValue *= ifelse(contains(stdevType, "REL"), value, SolarChem.unitconversionfactor(stdevUnit, io=devnull)) 
                stdevValue *= ifelse(contains(stdevType, "2"), 0.5, one(stdevValue)) 
                    # uncertaintynames = ("1S",  "REL", "2S", "2S-ABS", "S-REL")

                if variable ∈ keys(oxidelookup)
                    el = oxidelookup[variable]
                    oxc =  SolarChem.oxideconversion[el] 
                    el = ifelse(el==:Fe3, :Fe, el)
                    meas[el][ind] =  value * oxc
                    uncs[el][ind] = stdevValue * oxc
                elseif variable ∈ keys(ellookup)
                    el = ellookup[variable]
                    meas[el][ind] = value
                    uncs[el][ind] = stdevValue
                elseif variable ∉ ignored
                    push!(ignored,variable)
                end
            end
        end

        newrl = string(url,"&afterKey=", jsn[:afterKey][:batchId])
        jsn = JSON3.read( String( take!( Downloads.download(newrl, io) ) ) )
        Term.Progress.update!(job)
        Term.Progress.render(pbar)
    end
    Term.Progress.stop!(pbar)
    

    indout = Base.OneTo(ind)
    
    dataout = Vector{Vector{Float64}}(undef,2*length(meas))
    namesout = Vector{Symbol}(undef, length(dataout))

    ii=0
    for i = eachindex(meas)
        namesout[ii+1], namesout[ii+2] = (i, Symbol(:s,i))
        dataout[ii+1] = meas[i][indout]
        dataout[ii+=2] = uncs[i][indout]
    end


    namesout = (:name, :group, :type, :taxon, :citation, :dataset, :comment, namesout...)
    dataout = (name[indout], group[indout], type[indout], taxon[indout], citation[indout], dataset[indout], comment[indout], dataout...)
    out = NamedTuple{namesout}(dataout)

    if save
        Serialization.serialize(save,)
        print("\nData saved locally to:")
        println(Term.@cyan Term.@bold filesavelocation)
    else
        println(Term.@cyan Term.@bold "\nData was not saved locally.")
    end

    return showskipped ? (out, ignored) : out
end


"""

    unitconversionfactor(str; io=stdout)

Returns a unit conversion factor to scale measurements by to convert to mass fraction abundance based on the unit given in the heading `str`. By default, notifies user if a unit is ignored. 

includes: wt% (%, %m/m), mg/g, μg/g (ppm), ng/g (ppb), pg/g
excludes: at%, vol%, per mil, Ma, etc...


"""
function unitconversionfactor(s::AbstractString; io::IO=stdout)

    f = ifelse(occursin("wt%",s), 0.01, NaN)
    f = ifelse(occursin("%",s), 0.01, f) # assume -> wt%
    f = ifelse(occursin("%m/m",s), 0.01, f) # assume -> wt%
    
    f = ifelse(occursin("mg/g",s), 1e-3, f)

    f = ifelse(occursin("μg/g",s), 1e-6, f)
    f = ifelse(occursin("ppm",s), 1e-6, f)
    
    f = ifelse(occursin("ppb",s), 1e-9, f)
    f = ifelse(occursin("ng/g",s), 1e-9, f)
    
    f = ifelse(occursin("pg/g",s), 1e-12, f)

    isnan(f) && printstyled(io, "Caution: No concentration unit identified for `$s` [excluded from dataset, see ?unitconversionfactor for accepted units]\n", color=:yellow)

    f
end


######################
##### DEPRECATED #####
######################


"""

```julia
    SolarChem.loadastromatdata(file::String)
```

DEPRECATED: use [`astromatdata`](@ref), which directly accesses the AstroMat API.

Load data exported from Astromat as a csv. Given underlying group and type assignment functions, this will only load chondrite data. 

see also: [`assigngroup`](@ref), [`assigntype`](@ref), [`unitconversionfactor`](@ref), [`readvalue`](@ref)

"""
function loadastromatdata(file::String=string(@__DIR__,"/../data/astromat/astromat-6_27_2024-eventId-2294.csv"))

    x = DelimitedFiles.readdlm(file,',')
    # headings: "sample", "sample_url", "sample_url", "dataset_url", "dataset", "citation_url", "dataset_url", "citation", "citation_url", "collectionTypeName", "taxonType", "taxonName", "analysisType", "analyzedMaterialName", "analysisComment", "calcAvg", "numberOfReplicates"

    headers = string.(x[1,:])
    taxonName = x[:,findfirst(x -> x=="taxonName", headers)]
    # analyzedMaterialName = x[2:end,findfirst(x -> =="analyzedMaterialName", headers)] all "WHOLE ROCK"

    group_ = SolarChem.assigngroup.(taxonName)
    ni = findall(.!isempty.(group_))
    n = length(ni)

    group = group_[ni]
    type::Vector{Tuple{Vararg{Int}}} = SolarChem.assigntype.(taxonName[ni])
    name = string.(x[ni,findfirst(x -> x=="sample", headers)])
    comment = string.(x[ni, findfirst(x -> x=="analysisComment", headers)])
    citation = string.(x[ni,findfirst(x -> x=="citation_url", headers)])
    dataset = string.(x[ni,findfirst(x -> x=="dataset_url", headers)])

    out = (; name, group, type, comment, citation, dataset)

    oxides = (Ca = "CaO ", Na = "Na2O ", Ni = "NiO ", Co = "CoO ", Mg= "MgO ", Fe = "FeO ", Fe3 = "Fe2O3 ", Si = "SiO2 ", Al = "Al2O3 ", Mn = "MnO ", Cr = "Cr2O3 ", Ti ="TiO2 ", K = "K2O ", P = "P2O5 ")

    uncertaintynames = ("1S",  "REL", "2S", "2S-ABS", "S-REL")

    hbv = BitVector(undef,length(headers))

    for i = eachindex(SolarChem.periodictable())
        includesuncs=false # reset whether there are corresponding uncertainties.

        el = SolarChem.periodictable()[i]

        v = fill(NaN,n)
        sv = copy(v)
        elstr = SolarChem.strelements[i]

# Identify all measurement and uncertainty headers
        @inbounds for ii = eachindex(hbv)
            header = headers[ii]
            b  = contains(header, elstr)
            if el ∈ keys(oxides)
                b |= contains(header, oxides[el])
                if el == :Fe # include both iron oxides
                    b |= contains(header,oxides.Fe3)
                end
            end

            if el == :O # ensure oxides are not counted as O measurements
                @inbounds for ox in oxides
                    b &= !contains(header, ox)
                end
            end
# Exclude all elemental headers with lab, comment, or method descriptors --- these do not contain numerical chemical data.
            b &= !contains(header, "lab")
            b &= !contains(header,"method")
            b &= !contains(header,"comment")
            hbv[ii] = b
        end 

        heads = findall(hbv)

        @inbounds for j ∈ eachindex(heads)

            ij = heads[j]

            header = headers[ij]
            xx = view(x,:,ij)

            cf = SolarChem.unitconversionfactor(header)

# Convert from oxide to elemental.
            if el ∈ keys(oxides) && contains(header,"O")
                cf *= ifelse(contains(header, "Fe2O3"), oxideconversion.Fe3, oxideconversion[el])
            end

# If header is not an uncertainty, just read values as numbers or nans.
# else if header is an uncertainty, identify it and calculate it.
            if iszero(sum(contains.(header,uncertaintynames))) # if header is not an uncertainty class....

                @inbounds for k = 1:n
                    r = SolarChem.readvalue(xx[ni[k]]) * cf
                    vk = v[k]
                    if el == :Fe # Multiple measurements of Fe may need to sum
                        vk += ifelse(isnan(r),0,r)
                        vk = ifelse(iszero(vk),NaN,vk)
                    end
                    v[k] = ifelse(isnan(vk), r, vk) 
                end
            else 
                includesuncs = true
                uncfactor = ifelse(occursin("2S",header),0.5,1)
# For relative uncertainties
                if contains(header,"REL")
                    @inbounds for k = 1:n
                        r = SolarChem.readvalue(xx[ni[k]]) * v[k] * uncfactor # requires that v has already been filled, for this cell, which should be the case
                        svk = sv[k]  
                        sv[k] = ifelse(isnan(r), svk, r) # overwite whatever's there with the new relative uncertainty.
                    end
# For absolute uncertainties
                else
                    @inbounds for k = 1:n
                        r = SolarChem.readvalue(xx[ni[k]]) * cf * uncfactor
                        svk = el == :Fe ? sqrt(sv[k]^2 + ifelse(isnan(r),0,r*r)) : sv[k] # for conditions where there might already be a value there (Fe), add in quadrature
                        sv[k] = ifelse(isnan(svk), r, svk) 
                    end
                end
            end
        end
        # Add new data to output NamedTuple, if there is any:
        if countnotnans(v) > 0
            if includesuncs
                out = (; zip((keys(out)...,el, Symbol(:s,el)), (out..., v, sv))...)
            else
                out = (; zip((keys(out)..., el), (out..., v))...)
            end
        end
    end 
    return out
end



## Supporting functions for `loadastromatdata`

"""

    readvalue(x)

Return a value `x` that is either a `Number` as a float or parses `x` if it's a String. An empty String yields `NaN`. Tolerates usage of commas in numbers.

"""
readvalue(x::Number) = float(x)
readvalue(x::AbstractString) = parse(Float64,ifelse(isempty(x), "NaN", replace(x,',' => "")))




"""

    assigntype(str)

    Assign all corresponding genetic groups to a String `str` containing meteorite group/type information. See [`innergroups`](@ref) and [`outergroups`](@ref) for groups.

"""
function assigntype(v::AbstractString)
    
    i = findfirst(isdigit,v)

    if isnothing(i)
        i = 1
        v = ifelse(contains(v,"IMP"), "7", "0")
    end

    c = c2 = parse(Int,v[i])

    if contains(v, '-') 
        ic2 = findfirst(x -> x=='-',v )+1
        isdigit(v[ic2]) && (c2 = parse(Int,v[ic2]))
    end

    if contains(v,'/')
        ic2 = findlast(x -> x=='/',v)+1
        if isdigit(v[ic2])
            c2 = parse(Int,v[ic2])
        end
    end   
    (c:c2...,)
end



"""

```julia
assigngroup(str)
```

Assign all corresponding genetic groups to a String `str` containing meteorite group/type information. See [`innergroups`](@ref) and [`outergroups`](@ref) for groups.

"""
function assigngroup(v::AbstractString)
    s = ()
    s = ifelse(isEH(v), (s..., :EH), s)
    s = ifelse(isEL(v), (s..., :EL), s)
    s = ifelse(isE(v), (s..., :E), s)
    s = ifelse(isH(v), (s..., :H), s)
    s = ifelse(isL(v), (s..., :L), s)
    s = ifelse(isLL(v), (s..., :LL), s)
    s = ifelse(isR(v), (s..., :R), s)
    s = ifelse(isCI(v), (s..., :CI), s)
    s = ifelse(isCV(v), (s..., :CV), s)
    s = ifelse(isCM(v), (s..., :CM), s)
    s = ifelse(isCR(v), (s..., :CR), s)
    s = ifelse(isCH(v), (s..., :CH), s)
    s = ifelse(isCB(v), (s..., :CB), s)
    s = ifelse(isCK(v), (s..., :CK), s)
    s = ifelse(isCL(v), (s..., :CL), s)
    s = ifelse(isCO(v), (s..., :CO), s)
    s = ifelse(isCung(v), (s..., :Cung), s)
    s
end

## Chondrite group is_ functions:

isEH(v::AbstractString) = contains(v, "EH")

isEL(v::AbstractString) = 0 < sum(contains.(v,("EL3", "EL3.", "EL4", "EL5", "EL6", "EL7", "EL/", "/EL", "EL(", "(EL)"))) 

function isE(v::AbstractString)
    c = ("E3", "E3.", "E4", "E5", "E6", "E7", "ENST")
    b = 0 < sum(contains.(v, c)) 
    b &= !contains(v,"ACHON") # exclude aubrites
end

function isH(v::AbstractString)
    c = ("H3", "H3.", "H4", "H5", "H6", "H7", "H/", "/H", "H-")
    b = 0 < sum(contains.(v, c)) 
    b *= !contains(v,"EH") 
    b *= !contains(v,"CH") #none in astromat when built but this could always change
end

function isL(v::AbstractString)
    c = ("L3", "L3.", "L4", "L5", "L6", "L7", "L/", "/L", "L(", "(L)")
    b = 0 < sum(contains.(v, c))  
    b *= !contains(v,"EL") 
    b *= !contains(v,"CL") #none in astromat when built but this could always change
    b *= !contains(v,"ACAPULCOITE")
    b *= !contains(v,"LODRANITE")
    b *= ((contains(v,"/") | contains(v,")")) | !contains(v,"LL"))
end

isLL(v::AbstractString) = 0 < sum(contains.(v, ("LL3", "LL3.", "LL4", "LL5", "LL6", "LL7", "LL/", "/LL", "LL(", "(LL)")))

isR(v::AbstractString) = 0 < sum(contains.(v, ("R3", "R3.", "R4", "R5", "R6", "R7", "R/", "R(", "(R)")))

isCI(v::AbstractString) = contains(v, "CI1") 

isCV(v::AbstractString) = contains(v, "CV")

isCM(v::AbstractString) = contains(v, "CM") & !contains(v, "EUCRITE")

isCR(v::AbstractString) = contains(v, "CR") & !contains(v, "EUCRITE")

isCH(v::AbstractString) = 0 < sum(contains.(v, ("CH1", "CH2", "CH3")))

isCB(v::AbstractString) = contains(v, "CB")

isCK(v::AbstractString) = 0 < sum(contains.(v, ("CK3", "CK4", "CK5", "CK6", "CK7")))

isCL(v::AbstractString) = contains(v, "CL") & !contains(v, "UNCL")

isCO(v::AbstractString) =  0 < sum(contains.(v, ("CO1", "CO2", "CO3")))

isCung(v::AbstractString) = 0 < sum(contains.(v,("C-UNG", "C1", "C2", "C3", "C4")))



# unique(taxonName) =["LL3 CHONDRITE", "H4", "H3", "H6", "C-UNG", "H5", "L3", "CM", "DIOGENITE", "EUCRITE-MMICT", "EUCRITE", "H3.3 CHONDRITE", "CH3", "H3.0", "CK5", "EL3 CHONDRITE", "CM2", "ANGRITE", "EUCRITE-CM", "ACHONDRITE-UNG", "CR2", "", "UREILITE", "R4", "STONE-UNCL", "LUNAR", "EH4 CHONDRITE", "CV3", "CO3.5", "R3.8-5", "CO3.7", "L6", "E4", "EL5 CHONDRITE", "CI1", "L/LL4", "IRON-IVA", "LL5 CHONDRITE", "LL3.2/3.4 CHONDRITE", "EUCRITE (POLYMICT)", "H6 CHONDRITE", "L6 CHONDRITE", "L5 CHONDRITE", "CO3.6 CHONDRITE", "SHERGOTTITE", "L3.5 CHONDRITE", "L3.4 CHONDRITE", "L3.7-3.9 CHONDRITE", "MESOSIDERITE", "L3-6 CHONDRITE", "IRON-UNGROUPED", "EH3 CHONDRITE", "H3.7 CHONDRITE", "LL3.7 CHONDRITE", "CO3.0 CHONDRITE", "H3.9 CHONDRITE", "AUBRITE", "NIPR COLLECTION", "EUCRITE (UNBRECCIATED)", "CV3 CHONDRITE", "LUNAR-ANORTH. BRECCIA", "EL6 CHONDRITE", "ACAPULCOITE", "CM2 CHONDRITE", "CO3.4 CHONDRITE", "H5 CHONDRITE (IN ICE)", "H5 CHONDRITE", "UREILITE (AUG-BEARING)", "CK4 CHONDRITE", "L4 CHONDRITE", "LL3.2/3.5 CHONDRITE", "L3.4-3.7 CHONDRITE", "LL3.3 CHONDRITE", "LL6 CHONDRITE", "CM1/2 CHONDRITE", "CO3.5 CHONDRITE", "MARTIAN (OPX)", "SNC ORTHOPYROXENITE", "H4 CHONDRITE", "BRACHINITE", "L3.8 CHONDRITE", "LL3.4 CHONDRITE", "IRON-IIIAB", "L3.9 CHONDRITE", "EUCRITE (MG-RICH)", "H3.5-4 CHONDRITE", "CH3 CHONDRITE", "R3.6 CHONDRITE", "H4/5", "H3.5 CHONDRITE", "H3.4 CHONDRITE", "H5/6", "CHONDRITE-UNCL", "L3.7 CHONDRITE", "IRON, IAB-MG", "EUCRITE-PMICT", "L5", "CO3.6", "MARTIAN (SHERGOTTITE)", "L3.6 CHONDRITE", "MESOSIDERITE-B1", "L4", "CO3.0", "HOWARDITE", "L3.2 CHONDRITE", "LL3.5 CHONDRITE", "L3.7-4 CHONDRITE", "EUCRITE-UNBR", "CV3-AN", "LUNAR (ANORTH)", "H3.6 CHONDRITE", "LUNAR (GABBRO)", "EH4", "C2-UNG", "LL4 CHONDRITE", "EH4/5", "LL3.15", "LL3.9", "H/L3.9", "OC", "LL6 CHON. (BRECCIA)", "CR2 CHONDRITE", "CO3 CHONDRITE", "L/LL6", "CK3", "IRON-IVB", "CO3", "R3.8", "CK4", "MARTIAN (CHASSIGNITE)", "R3.5-6", "C3-UNG", "CK4/5", "EL4 CHONDRITE", "UREILITE-PMICT", "R3-4", "CK3-AN", "H-IMP MELT", "H3.8 CHONDRITE", "H~6", "LUNAR-BASALT", "H5-6 CHONDRITE", "CO3 CHONDRITE (ANOMALOUS)", "H5-7", "IRON, IIAB", "LUNAR (FELDSP. BRECCIA)", "UREILITE (POLYMICT)", "CK5 CHONDRITE", "IRON-IAB", "CM1 CHONDRITE", "C2 CHONDRITE UNGROUPED", "IRON-IIE (ANOMALOUS)", "L3.3-3.6 CHONDRITE", "ACAPULCOITE/LODRANITE", "IRON-IAB (ANOMALOUS)", "LUNAR-BASALTIC BRECCIA", "EUCRITE (BRECCIATED)", "UREILITE (POLYMICT ?)", "CK5/6 CHONDRITE", "L3.0 CHONDRITE", "CV3 CHONDRITE (REDUCED)", "H4 CHONDRITE (ANOMALOUS)", "EH4/5 CHONDRITE", "L/LL3.2 CHONDRITE", "TERRESTRIAL ROCK", "IRON-IIE", "CO3.3", "LODRANITE", "CR", "ACHON. UNGROUPED", "DIOGENITE (OLIVINE)", "EL4/5", "IRON-OCTAHEDRITE", "CHONDRITE UNGROUPED", "DIOGENITE (UNIQUE)", "CM CHONDRITE (ANOMALOUS)", "CR1 CHONDRITE", "LL3.8 CHONDRITE", "IRON, IAB-UNG", "L5-6", "C4-UNG", "WINONAITE", "CBB", "EL6/7", "L3.7-6", "R3-6", "EL7", "L/LL3.4", "CO3.8", "CO3.2", "K3", "L/LL5", "LL3.2", "R CHONDRITE", "H CHONDRITE (IMPACT MELT)", "L CHONDRITE (IMPACT MELT)", "R4 CHONDRITE", "CM1-2 CHONDRITE", "AUBRITE (ANOMALOUS)", "R6 CHONDRITE", "LUNAR (BASALT)", "LL5", "UREILITE (ANOMALOUS)", "CK6 CHONDRITE", "C3 CHONDRITE UNGROUPED", "L3.1 CHONDRITE", "L3.2-3.5 CHONDRITE", "E3 CHONDRITE (ANOMALOUS)", "CHON. (KAKANGARI-LIKE)", "EH5 CHONDRITE", "L5/6", "L(LL)3.05 CHONDRITE", "CV3 CHONDRITE (ANOMALOUS)", "R3 CHONDRITE", "L3.10 CHONDRITE", "NAKHLITE", "CB CHONDRITE", "LUNAR-FELDSPATHIC BRECCIA", "H5 CHONDRITE (ANOMALOUS)", "CR3 CHONDRITE", "CHONDRITE-UNG", "CM1", "CV2", "MARTIAN (NAKHLITE)", "IRON, IID-AN", "LL3.6 CHONDRITE", "CM-AN", "UNKNOWN", "MARTIAN (POLYMICT BRECCIA)", "R3/4", "R3", "CK5/6", "R3.9", "R5", "DIOGENITE-AN", "H", "ENST ACHON", "MARTIAN", "EUCRITE-AN", "DIOGENITE-PM", "LL5/6 CHONDRITE", "LL3.1", "LL3.05", "DIOGENITE-OLIVINE", "R3.7", "MARTIAN (AUGITE BASALT)", "H4-6 CHONDRITE", "C3", "CM1/2", "ACHONDRITE-PRIM", "E6", "CO3.4", "CH CHONDRITE", "CK3 CHONDRITE", "R3.8-6 CHONDRITE", "CBA", "ENSTATITE METEORITE UNGR", "CBB CHONDRITE", "CV3.4", "EUCRITE-BR", "H3.2-3.7", "MESOSIDERITE-AN", "L3.6-4 CHONDRITE", "R3.8-6", "R3-5", "KREEP BASALT", "H-METAL", "LL3.00", "L6/7", "H3-6", "C", "H/L3.6", "LL3.1-3.5 CHONDRITE", "H3-4", "LL7 CHONDRITE", "LODRANITE-AN", "H-AN", "LL3.0", "H7", "L4/5", "EH", "LL", "IRON, IIIAB", "IRON, IIE-AN", "H3/4", "E5-AN", "IRON, IAB COMPLEX", "EH6-AN", "H3.2-AN", "C1/2-UNG", "EH-IMP MELT", "CK6", "EH6 CHONDRITE", "PALLASITE, UNGROUPED", "R6"]




# from_astromat = [:Fe, :As, :Cu, :Mn, :Ca, :Mg, :Zn, :S, :Al, :FeOT, :NiO, :FeS, :Fe2O3, :SO2, :Cr, :Ti, :Na, :P, :CARBONATE, :Fe2P, :Fe3P, :CoO, :K, :Cd, :Fe2O3T, :LOI, :Si, :Ag, :MoO3, :SO3, :V2O3, :BaO, :ST, :TiO, :ZrO2, :Ir, :Ga, :Ni, :Au, :Cr, :Co, :As, :Ge, :W, :Re, :Pt, :Cu, :Sb, :Sc, :K, :Zn, :V, :Se, :Os, :Hf, :Na, :Sr, :Ca, :C, :N, :Fe, :Mn, :Ru, :Rb, :Th, :U, :Zr, :Br, :Al, :Mg, :Ta, :Cs, :Pd, :Ti, :Te, :Y, :Cd, :Nb, :Mo, :Li, :Ag, :Tl, :Bi, :In, :Pb, :Sn, :P, :I, :S, :Si, :Be, :Rh, :B, :Hg, :Ba, :Sm, :Lu, :Yb, :La, :Eu, :Nd, :Cd, :Tb, :Ba, :Dy, :Ho, :Gd, :Er, :Pr, :Tm]

# elements = [:Li, :Be, :B, :C, :N, :Na, :Mg, :Al, :Si, :P, :S, :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Rb, :Sr, :Y, :Zr, :Nb, :Mo, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn, :Sb, :Te, :I, :Cs, :Ba, :La, :Pr, :Nd, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, :Hf, :Ta, :W, :Re, :Os, :Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Th, :U]

# extras = [:FeOT, :FeS, :Fe2O3, :SO2, :CARBONATE, :Fe2P, :Fe3P, :Fe2O3T, :LOI,:MoO3, :SO3, :V2O3, :BaO, :ST, :TiO, :ZrO2, ]

# volatiles = [:H, :CO2, :F, :H2OM, :H2O, :H2OP, :N2, :O, :HCl, :O2, :CO1  ]
# volunits = ["umol/g", "ccstp", "ccstp/g"]