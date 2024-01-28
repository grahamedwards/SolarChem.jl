"""

    estimateuncertainty(d, unc; maxpctunc, minuncs)

Calculate an assumed uncertainty `unc` (in %) for all analytes in NamedTuple `d`. If there are `> minuncs` (30 by default) uncertainties corresponding to measurements of an analyte, it calculates unreported uncertainties as the mean relative uncertainty. Replaces any values exceeding the maximum percent uncertainty (`maxpctunc`) with `unc`. 

Requires that uncertainties are given as 1σ with keys `:sX` where X is the analyte name (e.g. `:Na` and `:sNa`).

"""
function estimateuncertainty(d::NamedTuple, unc::Number; maxpctunc::Number=50., minuncs::Int=30)
    unc = .01float(unc)
    maxunc = .01float(maxpctunc)
    @assert unc < maxunc

    extantuncs = keys(d)[findall(x -> startswith(string(x),"s"),keys(d))]
    analytes = keys(d)[findall(x -> !startswith(string(x),"s"),keys(d))]

    @inbounds for k in analytes
        sk = Symbol(:s,k)
        u = sk ∈ extantuncs ?  d[sk] : fill(NaN,length(d[k]))

        unc = count(!isnan,d[k]) > minuncs ? nanmean(d[sk]./d[k]) : unc

        @inbounds @simd for i = eachindex(d[k])
            x, sig = d[k][i], u[i]
            rsig = ifelse(isnan(sig), unc, sig/x ) # convert NaN to assumed uncertainty
            u[i] = x * ifelse( rsig > maxunc, unc, rsig)
        end
        if sk ∉ extantuncs 
            println(sk)
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



#####

# Add a function that removes vector indices based on certain conditions. 

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