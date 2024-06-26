## Resampling functions

"""

```julia
bsresample(n, x, σ; w, rng)
```

Returns a Vector of `n` random samples from dataset `x` with (normally distributed) 1σ uncertainties `σ`. Optionally provide weights `w` (obtained from [`calcweights`](@ref), unweighted by default) and a (pseudo)random number generator `rng` (default: Xoshiro256++).

"""
bsresample(n::Int, data::Vector{T},sigma::Vector{T}; w::Vector{T}=[],  rng::Random.AbstractRNG=Random.Xoshiro()) where T <: Float64 = bsresample!(Vector{eltype(data)}(undef,n), data, sigma, w=w, rng=rng)



"""

```julia
bsresample!(v, x, σ; w, rng)
```

In-place version of `bsresample` that takes a vector `v` to fill with resamples.

see also: [`bsmean`](@ref)

"""
function bsresample!(resampled::Vector{T}, data::Vector{T},sigma::Vector{T}; w::Vector{T}=[],  rng::Random.AbstractRNG=Random.Xoshiro()) where T<:Float64
        
    weights = isempty(w) ? ones(eltype(data),length(data)) : w

    @assert length(data) == length(sigma) == length(weights)
    
    wₛ = cumsum(weights)
    wₛ ./= last(wₛ)

    @inbounds @simd for i = eachindex(resampled)
        row = searchsortedfirst(wₛ, rand(rng))

        resampled[i] = data[row] + randn(rng)*sigma[row]
    end
    resampled 
end



"""

```julia
bsmean(n, x, σ; w, rng)
```

Returns a Vector of `n` means, each calculated from a random resampling (with replacement) of dataset `x` with (normally distributed) 1σ uncertainties `σ`. Optionally provide `weights` (obtained from [`calcweights`](@ref), unweighted by default) and a (pseudo)random number generator `rng` (default: Xoshiro256++).

"""
bsmean(n::Int, data::Vector{T},sigma::Vector{T}; w::Vector{T}=[],  rng::Random.AbstractRNG=Random.Xoshiro()) where T <: Float64 = bsmean!(Vector{eltype(data)}(undef,n), data, sigma, w=w, rng=rng)



"""

```julia
bsmean!(v, x, σ; w, rng)
```

In-place version of `bsmean` that takes a vector `v` to fill with resampled means.

see also: [`bsmean`](@ref)

"""
function bsmean!(means::Vector{T}, data::Vector{T},sigma::Vector{T}; w::Vector{T}=[],  rng::Random.AbstractRNG=Random.Xoshiro()) where T <: Float64
    
    weights = isempty(w) ? ones(eltype(data),length(data)) : w

    @assert length(data) == length(sigma) == length(weights)

    nd = length(data)
    wₛ = cumsum(weights) 
    wₛ ./= last(wₛ)

    @inbounds @simd for h = eachindex(means)
        μ = zero(eltype(means))
        for i = 1:nd
            row = searchsortedfirst(wₛ, rand(rng))
            μ += data[row] + randn(rng)*sigma[row]
        end
        means[h] = μ/nd
    end
    means
end

#   #   #   #   #   #   #   #   #   #   #
# Calculations in elemental batches:    #
#   #   #   #   #   #   #   #   #   #   #



"""

    bootstrapelements(n::Int, data::NamedTuple, elements; resamplemeans=true, weighted=true, rng)

Returns a NamedTuple of vectors of `n` bootstrap resampled data for each element in the Tuple `elements`. 
    
Resamples Monte Carlo'ed means by default. Declare `resamplemeans`=`false` to return resampled values. 
    
By default weights resampling by sample abundance (based on occurences of unique meteorite names (field `:name` in `data`). To remove weighting, declare `weighted`=`false`. 

see also: [`trimnans`](@ref), [`calcweights`](@ref), [`bsresample`](@ref), [`bsmean`](@ref)

"""
function bootstrapelements(n::Int, d::NamedTuple, els::Tuple{Vararg{Symbol}}; resamplemeans::Bool=true, weighted::Bool=true, rng::Random.AbstractRNG=Random.Xoshiro())
# Checks for informative errors
    k = keys(d)
    
    weighted && @assert :name ∈ k "data must contain a :name field if weighted=true"

    present = ()
    @inbounds for el in els
        if el ∈ k 
            present = (present...,el)
            sel = Symbol(:s,el)
            @assert  sel ∈ k "no corresponding uncertainty :$sel for :$el"
        else
            printstyled("Caution: input name :$el is not a name in the provided dataset. I am skipping this label.\n",color=:yellow)
        end
    end 

    els = present
    x = NamedTuple{els}(Vector{Float64}(undef,n) for i = eachindex(els))

    @batch for i in eachindex(els)
        
        el = els[i] 

        trimmed = trimnans(d,el, alsoinclude=ifelse(weighted,(:name,),()))
        weights = weighted ? calcweights(trimmed.name) : ones(length(trimmed[el]))

        if isempty(trimmed[el])
            @warn "No non-NaN ratios of :el."
            x[el] .= NaN
        elseif resamplemeans
            bsmean!(x[el], trimmed[el], trimmed[Symbol(:s,el)], w=weights, rng=rng )
        else
            bsresample!(x[el], trimmed[el], trimmed[Symbol(:s,el)], w=weights, rng=rng )
        end
    end
    x
end 



"""

    bootstrapratios(n::Int, data::NamedTuple, elements, divisor::Symbol; resamplemeans=true, weighted=true, rng)

Returns a NamedTuple of vectors of `n` bootstrap resampled data for each element in the Tuple `elements`, as a ratio of element/`divisor` (calculated prior to resampling).
    
Resamples Monte Carlo'ed means by default. Declare `resamplemeans`=`false` to return resampled values. 
    
By default weights resampling by sample abundance (based on occurences of unique meteorite names (field `:name` in `data`). To remove weighting, declare `weighted`=`false`. 

see also: [`trimnans`](@ref), [`calcweights`](@ref), [`bsresample`](@ref), [`bsmean`](@ref)

"""
function bootstrapratios(n::Int, d::NamedTuple, els::Tuple{Vararg{Symbol}}, divisor::Symbol;  resamplemeans::Bool=true, fractional::Bool=false, weighted::Bool=true, rng::Random.AbstractRNG=Random.Xoshiro())
    # Checks for informative errors
    k = keys(d)
    @assert divisor ∈ k "divisor $divisor is not in dataset"
    
    present = ()
    @inbounds for el in els
        if el ∈ k
            if el==divisor
                printstyled("Caution: :$el is the divisor of the ratio. I am excluding this from the provided numerators.\n",color=:yellow)
            else
                present = (present...,el)
                sel = Symbol(:s,el)
                @assert  sel ∈ k "no corresponding uncertainty :$sel for :$el"
            end
        else
            printstyled("Caution: input name :$el is not a name in the provided dataset. I am skipping this label.\n",color=:yellow)
        end
    end 
    els = present

    x = NamedTuple{els}(Vector{Float64}(undef,n) for i = eachindex(els))

    sdivisor = Symbol(:s,divisor)

    @batch for j in eachindex(els)

        el = els[j]
        sel = Symbol(:s,el)
        trimmed = trimnans(d,(el, divisor), alsoinclude=(:name,))
        weights = weighted ? calcweights(trimmed.name) : ones(length(trimmed.name))

        r = fractional ? 
            trimmed[el] ./ (trimmed[divisor] .+ trimmed[el]) : 
            trimmed[el] ./ trimmed[divisor] 

        sr = similar(r)

        tsd, td = trimmed[sdivisor], trimmed[divisor] # divisor
        tsn, tn = trimmed[sel], trimmed[el] # numerator
        
        @inbounds @simd ivdep for i in eachindex(sr)
            sxd, xd0 = tsd[i], td[i] # divisor, σ,μ
            sxn, xn = tsn[i], tn[i] # numerator σ,μ
            xd = ifelse(fractional, xd0+xn, xd0)
            sxd = ifelse(fractional, sqrt(sxd * sxd + sxn * sxn), sxd)
            xnf = sxn/xn
            xdf =  sxd/xd
            sr[i] = r[i]*sqrt(xnf*xnf + xdf*xdf)
        end 

        if isempty(trimmed[el])
            printstyled("Caution: No non-NaN ratios of :$el.",color=:yellow)
            x[el] .= NaN
        elseif resamplemeans
            bsmean!(x[el], r, sr, w=weights, rng=rng )
        else
            bsresample!(x[el], r, sr, w=weights, rng=rng )
        end
        fractional && fraction2ratio!(x[el]) 
    end
    x
end 