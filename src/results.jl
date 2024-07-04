"""

    ratiosummary(resampled_measurements, resampled_ratios; ratiocounts, minratios)

Calculate ratio mean and standard deviation (1σ) from `resampled_measurements` calculated by [`bootstrapelements`](@ref) and `resampled_ratios` calculated by[`bootstrapratios`](@ref). 

If no ratios were calculated for a field in `resampled_ratios` (all values `NaN`), μ/σ are calculated from the corresponding `resampled_measurements`.

Optionally, provide `ratiocounts` calculated by [`countratios`](@ref) and a minimum number of ratios `minratios` for each element required to use `resampled_ratios`. If the ratio counts for an element is less than `minratios`, μ/σ are calculated from the corressponding elements in `resampled_measurements`.

---

`NamedTuple(x)` for some `x = ratiosummary...` will convert the result to a NamedTuple of [`Composiiton`](@ref)s.

"""
function ratiosummary(rsmeas::NamedTuple, rsratio::NamedTuple; ratiocounts::NamedTuple=(;), minratios::Int=-1)

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
    (numerators,ms)
end

NamedTuple(x::Tuple{Vector{Symbol}, Matrix{Float64}}) = (; zip(x[1], [Composition(x[2][i,:]...) for i in axes(x[2],1)])...)