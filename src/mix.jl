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
normpdf(x::T, m::T, s::T) where T<:Float64 = (0.3989422804014327/s) * exp(-(x-m)*(x-m)/(2.0*s*s)) # 1/sqrt(2π) = 0.3989422804014327 



"""

    fastll(x, m, s)

Calculate the relative log-likelihood that `x` is drawn from the Normal distribution `m` ± `s`. `fastll` excludes the additive constant `-log(sqrt(2π)*σ)`, which can be multiplied to summations of `exp(fastll(..))` later for accurate likelihood calculations. All inputs must be `Float64`s.

"""
fastll(x::T, m::T, s::T) where T<:Float64 = -(x-m)*(x-m)/(2.0*s*s)


"""

    lldist(μ, σ, μₒ, σₒ; n=100)

Calculates the log-likelihood that the normally distributed measurement(s) `μ` ± `σ` was/were drawn from the model composition `μₒ` ± `σₒ`, integrated over `n` nodes of the 3σ limits of the two distributions. 

Accepts discrete `Float64` values or vectors of `Float64`s for `μ` and `σ`. Values of `NaN` are ignored and do not contribute to log-likelihood calculations.

see also: [`normpdf`](@ref)

"""
function lldist(m::Float64, s::Float64, mo::Float64, so::Float64; n::Int=100)   
    l = 0.
    nsigs = 3.
    siglims = m-nsigs*s, m+nsigs*s, mo-nsigs*so, mo+nsigs*so
    xmin, xmax = extrema(siglims)
    x = LinRange(xmin, xmax, n)
    @inbounds @simd for i = eachindex(x)
        l += @inline exp(fastll(x[i],m,s) + fastll(x[i], mo, so))
    end
    log(l * 0.3989422804014327 * (xmax-xmin)/(s * so * (n-1))) # scale by Δx [Δx = (xmax-xmin)/(n-1) ] and divide by sqrt(2π)*σ*σₒ
end

function lldist(m::Vector{Float64}, s::Vector{Float64}, mo::Float64, so::Float64; n::Int=100)
    @assert length(m) === length(s)
    ll = 0.
    if !isnan(mo) & !isnan(so)
        @inbounds for i= eachindex(m)
            mi, si = m[i], s[i]
            if !isnan(mi) & !isnan(si)
                ll += lldist(mi, si, mo, so,n=n)
            end
        end
    end
    ll
end



"""

    jump(p::Fractions, j::Fractions; rng)

Perturb a random field of p (`outer`,`sun`) with a random gaussian jump with σ described by the corresponding field in `j`. Ensures that the new proposal satisfies that `0 < outer < 1` and `0 < sun < 1`. Optionally supply a random number generation seed `rng`.

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

Calculate the solar photosphere-normalized composition of a model Sun that ingests `f.sun` M⊙ of chondritic material comprised of `f.outer` proportion of `outer` solar system composition and `1 - f.outer` proportion of `inner` solar system composition. Both `inner` and `outer` are vectors containing the μ and σ in the first and second inidices, respectively.

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
function solarmixmetropolis(chainsteps::Int, proposal::Fractions, jumpsize::Fractions, innersolarsystem::NamedTuple, outersolarsystem::NamedTuple; burnin::Int=0, solarunc::Bool=true, stars=solartwins())

    jumpscale=2.9

    iss = innersolarsystem
    oss = outersolarsystem

    @assert keys(iss) == keys(oss)
    ksol = keys(iss)
    kstars = keys(stars)
    add_el, pt = (), periodictable()
    @inbounds for i = kstars
        if i ∈ pt && i ∉ ksol
            add_el = (add_el...,i)
        end
    end
    addnt = NamedTuple{add_el}((Composition(NaN,NaN) for i = length(add_el)))
    iss, oss = (; iss..., addnt...), (; oss..., addnt...)
    ksol = keys(iss)

    s = sun(:Fe, logspace=false)

    p, j = proposal, jumpsize
    ϕ = p 
    llϕ = zero(Float64)
    @inbounds @batch reduction = ((+,llϕ)) for j = eachindex(ksol) # 
        el = ksol[j]
        slm = solarlogmix(iss[el], oss[el], s[el], ϕ, solarunc=solarunc)
        llϕ += lldist(stars[el],stars[Symbol(:s,el)], slm...)
    end
    ll=llϕ

    @inbounds for i = 1:burnin
        ϕ,jumpname,jumpval = jump(p,j;rng=rng)

        llϕ = zero(ll)
        @inbounds @batch reduction = ((+,llϕ)) for j = eachindex(ksol) # 
            el = ksol[j]
            slm = solarlogmix(iss[el], oss[el], s[el], ϕ, solarunc=solarunc)
            llϕ += lldist(stars[el],stars[Symbol(:s,el)], slm...)
        end

        if log(rand(rng)) < (llϕ-ll) 
            jumpup = jumpval*jumpscale
            j = Fractions(ifelse(jumpname==:outer,(jumpup,f.sun),(f.outer,jumpup))...)
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
        @inbounds @batch reduction = ((+,llϕ)) for j = eachindex(ksol) # 
            el = ksol[j]
            slm = solarlogmix(iss[el], oss[el], s[el], ϕ, solarunc=solarunc)
            llϕ += lldist(stars[el],stars[Symbol(:s,el)], slm...)
        end

        if log(rand(rng)) < (llϕ-ll) 
            n_acceptance += 1
            jumpup = jumpval*jumpscale
            j = Fractions(ifelse(jumpname==:outer,(jumpup,f.sun),(f.outer,jumpup))...)
            p, ll = ϕ, llϕ  # update proposal and log-likelihood
        end
        vouter[i], vsun[i], vll[i] = p.outer, p.sun, ll
    end
    (; outer=vouter, sun=vsun, ll=vll, acceptance=n_acceptance/chainsteps)
end