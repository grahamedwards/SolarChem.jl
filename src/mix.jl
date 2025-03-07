## Mixing model

## Requisite: Compositional Models of inner and outer (which I am close to!)

## 2 parameters:
    # Proportion of inner vs outer solar system (outer proportion seems reasonable).
    # Fraction of solar mass

## Model structure:
    # Calculate mixed composition by adding and propagating uncs in quadrature
    # Add to composition of Sun, scaled by fraction of solar mass

## MCMC
    # loglikelihood calculation of each element vs. every solar twin. 

"""

    normpdf(x, m, s)

Calculate the probability density of `x` given a normal distribution with mean `m` and standard deviation `s`. All inputs must be `Float64`s.

``\\mathcal{L}(\\mu \\pm \\sigma | x ) = f(x | \\mu \\pm \\sigma) = \\frac{e^{-\\frac{1}{2} \\left( \\frac{x-\\mu}{\\sigma}\\right)}}{\\sigma\\sqrt{2\\pi}}``

"""
normpdf(x::T, m::T, s::T) where T<:Float64 = (0.3989422804014327/s) * exp(fastll(x,m,s)) # 1/sqrt(2π) = 0.3989422804014327 



"""

    fastll(x, m, s)

Calculate the relative log-likelihood that `x` is drawn from the Normal distribution `m` ± `s`. `fastll` excludes the additive constant `-log(sqrt(2π)*σ)`, which can be multiplied to summations of `exp(fastll(..))` later for accurate likelihood calculations. All inputs must be `Float64`s.

"""
fastll(x::T, m::T, s::T) where T<:Float64 = -(x-m)*(x-m)/(2.0*s*s)


"""

    lldist(μ, σ, μₒ, σₒ; n=100)

Calculates the log-likelihood that the normally distributed measurement(s) `μ` ± `σ` was/were drawn from the model composition `μₒ` ± `σₒ`, integrated over `n` nodes of the 3σ limits of the model distribution. 

Accepts discrete `Float64` values or vectors of `Float64`s for `μ` and `σ`. Values of `NaN` in modeled data are ignored and do not contribute to log-likelihood calculations.

see also: [`normpdf`](@ref)

"""
function lldist(m::Float64, s::Float64, mo::Float64, so::Float64; n::Int=100)   
        likelihood = 0.
        nsigs = 3.0 * s 
        x = LinRange(m-nsigs, m+nsigs, n)
        Δx = step(x)
        @inbounds for xi in x
            likelihood += exp(SolarChem.fastll(xi, mo, so) + SolarChem.fastll(xi, m, s))
        end
        L = Δx * likelihood * 0.15915494309189535 /(s*so)
        log(ifelse(iszero(L), eps(), L)) # scale by Δx and divide by 2π*σ*σₒ (2π = sqrt(2π)^2)
end

function lldist(m::Vector{Float64}, s::Vector{Float64}, mo::Float64, so::Float64; n::Int=100)
    @assert length(m) === length(s)
    ll = 0.
    if !isnan(mo+so)
        @inbounds for i= eachindex(m)
            ll += lldist(m[i], s[i], mo, so,n=n)
        end
    end
    ll
end



"""

    jump(p::Fractions, j::Fractions; rng)

Perturb a random field of p (`outer`,`sun`) with a random gaussian jump σ described by the corresponding field in `j`. Ensures that the new proposal satisfies that `0 < outer < 1` and `0 < sun < 1`. Optionally supply a random number generation seed `rng`.

see also: [`Fractions`](@ref)

"""
function jump(p::Fractions, j::Fractions; rng::Random.AbstractRNG=Random.Xoshiro())
    jn = rand(rng,fieldnames(Fractions))
    jv = getfield(j,jn) * randn(rng)
    jp = getfield(p,jn) + jv
    
    while (0 > jp) | (jp > 1)
        # println(jp) # for testing
        jv = getfield(j,jn) * randn(rng)
        jp = getfield(p,jn) + jv
    end
    Fractions(ifelse(jn ==:outer, (jp, p.sun),(p.outer, jp))...), jn, jv
end



"""

    solarlogmix(inner, outer, solar, f::Fractions; solarunc=true)

Calculate the solar photosphere-normalized composition of a model Sun that ingests `f.sun` M⊙ of chondritic material comprised of `f.outer` proportion of `outer` solar system composition and `1 - f.outer` proportion of `inner` solar system composition. The args `inner`, `outer`, and `solar` are all instances of a [`Composition'](@ref) struct. 

The result is reported in dex, as the ratio of log-ratios: `log10(X/D) - log10(X⊙/D⊙)` where `D` denotes the ratio divisor (e.g. Fe, H) and ⊙ denotes solar compositions.

"""
function solarlogmix(i::T, o::T, solar::T, f::Fractions; solarunc::Bool=true) where T <: Composition

    finner = 1 - f.outer
    iomix = i.m*finner + o.m * f.outer
    siomix2 = i.s*i.s*finner*finner + o.s*o.s*f.outer*f.outer

    fssol = ifelse(solarunc, solar.m/solar.s, 0.)

    added = f.sun * iomix/solar.m
    sig = added * sqrt(fssol*fssol + siomix2 / (iomix*iomix))

    solmixed, sig = (added + 1.0, sig) ./ (1.0 + f.sun)

    log10(solmixed), 0.43429448190325176 * sig / solmixed # 0.43429448190325176 = 1/log(10)
end 


"""

    solarmixmetropolis...

Metropolis algorithm to estimate the requisite mixture of chondritic components added to the solar photosphere to match it to solar twin compositions.

"""
function solarmixmetropolis(chainsteps::Int, proposal::Fractions, jumpsize::Fractions, innersolarsystem::NamedTuple, outersolarsystem::NamedTuple; burnin::Int=0, solarunc::Bool=true, stars=solartwins(), rng::Random.AbstractRNG=Random.Xoshiro())

    jumpscale=2.9

    iss = innersolarsystem
    oss = outersolarsystem

# Ensure that any elements from kstars missing in ksol are added to i/oss as NaN compositions. This helps the markov chain stay flexible and run smoothly. Also flags NaN and 0 values in stars. 
    @assert keys(iss) == keys(oss)
    ksol = keys(iss)
    kstars = keys(stars)
    add_el, pt = (), periodictable()
    @inbounds for i = kstars
        if i ∈ pt 
            @inbounds for ii = eachindex(stars[i])
                iszero(stars[Symbol(:s,i)][ii]) && @warn "σ=0 value for $i (row $ii) in prior dataset. Unintended results likely."

                isnan(stars[Symbol(:s,i)][ii]) && @warn "σ=0 value for $i (row $ii) in prior dataset. Unintended results likely."

                isnan(stars[i][ii]) && @warn "NaN value for $i (row $ii) in prior dataset. Unintended results likely."

            end
    
            if i ∉ ksol
                add_el = (add_el...,i)
            end
        end
    end
    addnt = NamedTuple{add_el}((Composition(NaN,NaN) for i = eachindex(add_el)))
    iss, oss = (; iss..., addnt...), (; oss..., addnt...)
    ksol = keys(iss)

    s = sun(:Fe, logspace=false)

    p, j = proposal, jumpsize
    ϕ = p 
    llϕ = zero(Float64)
    @inbounds @batch reduction = ((+,llϕ)) for ii = eachindex(ksol)  
        el = ksol[ii]
        slm = solarlogmix(iss[el], oss[el], s[el], ϕ, solarunc=solarunc)
        llϕ += lldist(stars[el],stars[Symbol(:s,el)], slm...)
    end
    ll=llϕ

    @inbounds for i = 1:burnin
        ϕ,jumpname,jumpval = jump(p,j;rng=rng)

        llϕ = zero(ll)
        @inbounds @batch reduction = ((+,llϕ)) for ii = eachindex(ksol) 
            el = ksol[ii]
            slm = solarlogmix(iss[el], oss[el], s[el], ϕ, solarunc=solarunc)
            llϕ += lldist(stars[el],stars[Symbol(:s,el)], slm...)
        end

        if log(rand(rng)) < (llϕ-ll) 
            jumpup = abs(jumpval)*jumpscale
            j = Fractions(ifelse(jumpname==:outer,(jumpup,j.sun),(j.outer,jumpup))...)
            p, ll = ϕ, llϕ  # update proposal and log-likelihood  
        end
    end 

    n_acceptance = 0
    vll = Vector{Float64}(undef,chainsteps)
    vouter = similar(vll)
    vsun = similar(vll)

    @inbounds for i = 1:chainsteps
        ϕ,jumpname,jumpval = jump(p,j;rng=rng)
        llϕ = zero(ll)
        @inbounds @batch reduction = ((+,llϕ)) for ii = eachindex(ksol)
            el = ksol[ii]
            slm = solarlogmix(iss[el], oss[el], s[el], ϕ, solarunc=solarunc)
            llϕ += lldist(stars[el],stars[Symbol(:s,el)], slm...)
        end

        if log(rand(rng)) < (llϕ-ll) 
            n_acceptance += 1
            jumpup = abs(jumpval)*jumpscale
            j = Fractions(ifelse(jumpname==:outer,(jumpup,j.sun),(j.outer,jumpup))...)
            p, ll = ϕ, llϕ  # update proposal and log-likelihood
        end
        vouter[i], vsun[i], vll[i] = p.outer, p.sun, ll
    end
    (; outer=vouter, sun=vsun, ll=vll, acceptance=n_acceptance/chainsteps, lastjump=j)
end