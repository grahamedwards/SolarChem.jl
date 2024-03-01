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

    @inbounds for h = eachindex(means)
        μ = zero(eltype(means))
        @inbounds @simd for i = 1:nd
            row = searchsortedfirst(wₛ, rand(rng))
            μ += data[row] + randn(rng)*sigma[row]
        end
        means[h] = μ/nd
    end
    means
end


