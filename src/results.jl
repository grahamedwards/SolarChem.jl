"""

    ratiosummary(resampled_measurements, resampled_ratios; ratiocounts, minratios, NTout=false)

Calculate ratio mean and standard deviation (1σ) from `resampled_measurements` calculated by [`bootstrapelements`](@ref) and `resampled_ratios` calculated by[`bootstrapratios`](@ref). 

If no ratios were calculated for a field in `resampled_ratios` (all values `NaN`), μ/σ are calculated from the corresponding `resampled_measurements`.

Optionally, provide `ratiocounts` calculated by [`countratios`](@ref) and a minimum number of ratios `minratios` for each element required to use `resampled_ratios`. If the ratio counts for an element is less than `minratios`, μ/σ are calculated from the corressponding elements in `resampled_measurements`.

If `NTout=true`, returns a NamedTuple of [`Composition`](@ref)s. 

"""
function ratiosummary(rsmeas::NamedTuple, rsratio::NamedTuple; ratiocounts::NamedTuple=(;), minratios::Int=-1, NTout::Bool=false)

    divisor = rsratio.divisor

    @assert divisor ∈ keys(rsmeas)
    @assert length(rsmeas) === length(rsratio)


    numerators = collect(keys(rsratio))
    numerators = numerators[numerators .!= :divisor]

    rc = isempty(ratiocounts) ? NamedTuple{keys(rsratio)}((zeros(Int,length(rsratio))...,)) : ratiocounts
    
    dm, ds = Statistics.mean(rsmeas[divisor]), Statistics.std(rsmeas[divisor])
    dsdm2 = (ds/dm)^2


    ms = Matrix{eltype(rsmeas[1])}(undef,length(numerators),2)

    @inbounds for i = eachindex(numerators)
        x = numerators[i]
        @assert x ∈ keys(rsmeas) 

        ms[i,:] .= 
        if isnan(rsratio[x]) || rc[x] < minratios
            nm,ns = Statistics.mean(rsmeas[x]),Statistics.std(rsmeas[x])
            r = nm/dm 
            r, r*sqrt((ns/nm)^2 + dsdm2)
        else
            Statistics.mean(rsratio[x]), Statistics.std(rsratio[x])
        end
    end
    return NTout ? 
        (; zip(numerators, [Composition(ms[i,:]...) for i in axes(ms,1)])...) : 
        (numerators,ms)
end



"""

    removefrom(s, x)

Remove key `s` and corresponding value from NamedTuple `x`. Alerts user if `s` is not present in `x`.

---
```julia
julia> removefrom(:c, (; a=1, b=2, c=3))
(a = 1, b = 2)

julia> removefrom(:c, (; a=1, b=2))
`:c` is not present in this NamedTuple
(a = 1, b = 2)
```

"""
function removefrom(s::Symbol, x::NamedTuple)
    k = keys(x)
    if s in k
        kout = out = ()
        @inbounds for i in k
            if s != i
                kout = (kout..., i )
                out = (out..., x[i])
            end
        end
        y = (; zip(kout, out)...)
    else
        println("`:$s` is not present in this NamedTuple")
        y = x
    end
    y
end