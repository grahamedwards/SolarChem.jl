"""

    estimateuncertainty(d, unc; uncextrema, minuncs)

Calculate an assumed uncertainty `unc` (in %) for all analytes in NamedTuple `d`. Replaces any values outside the bounds of the minimum and maximum percent uncertainty (`uncextrema`= (0.01%, 50%) by default) with `unc`. If there are `> minuncs` (30 by default) uncertainties corresponding to measurements of an analyte, it calculates unreported uncertainties as the mean relative uncertainty of all measurements within `uncextrema`.

Requires that uncertainties are given as 1σ with keys `:sX` where X is the analyte name (e.g. `:Na` -> `:sNa`).

Using `estimateuncertainty!` overwrites `d`.

"""
estimateuncertainty(d::NamedTuple, unc::Number; uncextrema::Tuple{Number,Number}=(.1, 50.), minuncs::Int=30) = estimateuncertainty!(deepcopy(d), unc, uncextrema=uncextrema, minuncs=minuncs)

function estimateuncertainty!(d::NamedTuple, unc::Number; uncextrema::Tuple{Number,Number}=(.1, 50.), minuncs::Int=30)
    @assert uncextrema[1] < unc < uncextrema[2] "Ensure assumed uncertainty is within defined range of uncextrema"
    
    unc = 0.01float(unc)
    minunc, maxunc = 0.01 .* float.(uncextrema)

    extantuncs = keys(d)[findall(x -> startswith(string(x),"s"),keys(d))]
    analytes = keys(d)[findall(x -> !startswith(string(x),"s") && x ∈ SolarChem.periodictable(),keys(d))]

    @inbounds for k in analytes
        sk = Symbol(:s,k)
        u = sk ∈ extantuncs ?  d[sk] : fill(NaN,length(d[k]))

        unc = if count(!isnan,u) > minuncs 
            rsigs = d[sk]./d[k]
            vmean(rsigs[minunc .< rsigs .< maxunc])
        else 
            unc 
        end

        @assert unc < maxunc "Mean uncertainty ($(100*unc)%) exceeds maximum allowed ($(uncextrema[2])%). Inspect $sk data. "
        @assert unc > minunc "Mean uncertainty ($(100*unc)%) less than minimum allowed ($(uncextrema[1])%). Inspect $sk data. "

        @inbounds @simd for i = eachindex(d[k])
            x, sig = d[k][i], u[i]
            rsig = ifelse(isnan(sig), unc, sig/x ) # convert NaN to assumed uncertainty
            u[i] = x * ifelse(minunc < rsig < maxunc, rsig, unc)
        end
        if sk ∉ extantuncs 
            d = (; zip((keys(d)...,sk), (d..., u))...) 
        end
    end
    d
end
    
#= MetBase: minimum number of reported uncertainties to calculate an average from:
30 --> ["Li", "Na", "Mg", "Al", "S", "Ca", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Zn", "Yb", "Os"]
    All reasonable and σ<10%, except Yb and Os at σ~16%, but these don't matter in stellar ?

    Li, Cr, Zn all have >40 recorded uncs.
=#



"""

    trimnans(data, names; alsoinclude = SolarChem.metadata)

Given a NamedTuple of `data`, return only the values and uncertainties of all non-NaN instances of `name` --- a name (::Symbol) or `Tuple` of names from `data`. Also returns the the Vectors corresponding to the fields listed in `alsoinclude`, which is by definition those defined by [`SolarChem.metadata`](@ref).


"""
function trimnans(d::NamedTuple, kz::NTuple{Nkz,Symbol}; alsoinclude::NTuple{Nai,Symbol}=SolarChem.metadata()) where {Nkz, Nai}

    @inbounds for k in kz # make sure that uncertainties are included in `outkeys`
        @assert k ∈ keys(d) "$k is not a name in the supplied dataset."
        sk = Symbol(:s,k)
        sk ∈ keys(d) && (kz = (kz...,sk))
    end
    
    n = length(d[kz[1]])
    notnans = trues(n)
    
    @inbounds for k in kz
        @assert length(d[k]) == n "length of $k is not consistent"
        @inbounds @simd for i=1:n
            notnans[i] *= ifelse(isnan(d[k][i]), false, true) 
        end
    end
    outkeys = (alsoinclude..., kz...)
    (; zip(outkeys, (d[k][notnans] for k in outkeys))...)
end

trimnans(d::NamedTuple, k::Symbol; alsoinclude::Tuple=SolarChem.metadata()) = trimnans(d, (k,), alsoinclude=alsoinclude)



"""

    trimextremes(data; min=0, max=1)

Cycle through all elements in `data` and replace any value in excess of `min` and `max` (set by default to the physical limits of fractional mass abundance) to a NaN.

Using `trimextremes!` overwites `data`. 

"""
trimextremes(d::NamedTuple; min::Number=0., max::Number=1.) = trimextremes!(deepcopy(d),min=min,max=max)

function trimextremes!(d::NamedTuple; min::Number=0., max::Number=1.)
   
    counter = 0

    @inbounds for el in periodictable()
        if el ∈ keys(d)
            x = d[el]
            @inbounds @simd for i = eachindex(x)
                xi = x[i]
                keep = (min <= xi <= max)
                counter += ifelse(isnan(xi),0,!keep)
                x[i] = ifelse(keep, xi, NaN)
            end
        end
    end
    
    println("\n$counter extremes replaced with NaNs\n")
    return d
end



"""

    pulltopic(data, name, topic[s]; exactmatch=true)

Given a NamedTuple of `data`, returns a NamedTuple with all the names of `data` and only including those rows within under `name` (::Symbol) that contain the provided `topic`(`s`) --- as a String or Tuple of Strings.

"""
function pulltopic(d::NamedTuple, k::Symbol, topics::NTuple{Nt,String}; exactmatch::Bool=false) where {Nt}
    @assert k ∈ keys(d) "$k is not a name in the supplied dataset."
    
    n = length(d[k])
    keepers = falses(n)
    @inbounds for t=1:Nt
        @inbounds @simd for i=1:n
            x= d[k][i]
            y= topics[t]
            keepers[i] |= ifelse(exactmatch,isequal(x,y),contains(x,y))
        end
    end
    
    iszero(keepers) && @warn "No matches for given topic(s)"
    
    (; zip(keys(d), (d[k][keepers] for k in keys(d)))...)
end
pulltopic(d::NamedTuple, k::Symbol, c::String; exactmatch::Bool=false) = pulltopic(d,k,(c,), exactmatch=exactmatch)



"""

    pullgroup(data, group; exactmatch=false)

Given a NamedTuple of chondrite `data`, returns a NamedTuple with all the names of `data`, including only the rows corresponding to the provided `group`(s) --- as a Symbol or Tuple of Symbols.

"""
function pullgroup(d::NamedTuple, topics::NTuple{N,Symbol}; exactmatch::Bool=false) where {N}

    keepers = falses(length(d.group))
    @inbounds @simd for t = 1:N
        @inbounds @simd for i = eachindex(keepers)
            x= d.group[i]
            keepers[i] |= ifelse(exactmatch, topics == x , topics[t] ∈ x )
        end
    end
    iszero(keepers) && @warn "No matches for given group(s)"
    (; zip(keys(d), (d[k][keepers] for k in keys(d)))...)
end
pullgroup(d::NamedTuple, c::Symbol; exactmatch::Bool=false) = pullgroup(d, (c,), exactmatch=exactmatch)


"""

    pulltype(data, type; exactmatch=false)

Given a NamedTuple of chondrite `data`, returns a NamedTuple with all the names of `data`, including only the rows corresponding to the provided `type`(s) --- as an Integer or Tuple of Integers.

"""
function pulltype(d::NamedTuple, topics::NTuple{N,Int}; exactmatch::Bool=false) where {N}
    
    n = length(d.type)
    keepers = falses(n)
    @inbounds @simd for t = 1:N
        @inbounds @simd for i = 1:n
            x = d.type[i]
            keepers[i] |= ifelse(exactmatch, topics == x , topics[t] ∈ x )
        end
    end
    iszero(keepers) && @warn "No matches for given type(s)"
    (; zip(keys(d), (d[k][keepers] for k in keys(d)))...)
end
pulltype(d::NamedTuple, c::Int; exactmatch::Bool=false) = pulltype(d, (c,), exactmatch=exactmatch)



"""

    exclude(data::NamedTuple, name::Symbol, s)

Exclude all indices containing string (or Tuple of strings) `s` from the specified `name` in `data`.

"""
function exclude(d::NamedTuple, k::Symbol, s::NTuple{N, String}) where N
    @assert k ∈ keys(d) "$k is not a name in the supplied dataset."

    keep = trues(length(d[k]))
    x = d[k]
    @inbounds for i = 1:N 
        si = s[i]
        @inbounds @simd for i = eachindex(keep)
            keep[i] &= !contains(x[i],si)
        end
    end
    return (; zip(keys(d), (d[k][keep] for k in keys(d)))...)
end
exclude(d::NamedTuple, k::Symbol, s::String) = exclude(d,k,(s,))

#####

"""

    calcweights(v ; weightings)

Calculate weightings for a population `v` for a non-biased sampling of the population.  Optionally give `weightings` as a `Dict` of `"type" => x` to further suppress values of `"type"` by a factor of `x`.

e.g.
    julia> v = ["a", "b", "a"]
    julia> calcweights(v)
    3-element Vector{Float64}:
     0.5
     1.0
     0.5

"""
function calcweights(v::Vector; weights::Dict{String,T}=Dict{String,Float64}()) where T <: Number
    
    w  = Vector{Float64}(undef, length(v))
    iᵢ = BitVector(undef, length(v))
    uv = unique(v)
    uw = ones(length(uv))

    if length(weights)>0
        uk = unique(keys(weights))
        @inbounds @simd for i in eachindex(uv)
            ii = uv[i]
            ii ∉ uk && error("Missing key '$ii' in weights")
            uw[i] = 1/weights[ii]
        end
    end

    @inbounds for i in eachindex(uv)
        @inbounds @simd for j in eachindex(v)
            iᵢ[j] = (v[j] == uv[i]) # find all indices of in v
        end
        nᵢ = count(iᵢ) # Calculate counts of iᵗʰ sample
        x=  uw[i]/ nᵢ # Calculate weight from counts
        @inbounds @simd for j in eachindex(w)
            y = w[j]
            w[j] = ifelse(iᵢ[j], x,y) # find all indices of in v
        end
    end
    w
end