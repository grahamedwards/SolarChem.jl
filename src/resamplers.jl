## Resampling functions

## Bootstrap resampler function
# This assumes a gaussian analytical distribution

function bootstrapresample(n::Int64,data::Vector,sigma::Vector,weights::Vector; )
    @assert length(data) == length(sigma) == length(weights)
    wₛ = cumsum(weights)
    
    resampled = similar(data, n)
    #rows = Vector{Int}(undef, n)
    @inbounds @simd for i = 1:n
        row = searchsortedfirst(wₛ, rand())
        #rows[i] = row
        resampled[i] = data[row] + randn()*sigma[row]
    end
    resampled # (resampled,rows)
end

# Resample the dataset (randomly) with replacement
    # Calculate mean of new dataset, repeat "nrs" times
    # This assumes a gaussian analytical distribution
function bootstrapmean(n::Int64,data::Vector,sigma::Vector,weights::Vector)
    @assert length(data) == length(sigma) == length(weights)
    nd = length(data)
    wₛ = cumsum(weights) 
    wₛ ./= last(wₛ)
    means = Vector{eltype(data)}(undef,nrs)
    @inbounds for h = 1:n
        μ = zero(eltype(means))
        @inbounds for i = 1:nd
            row = searchsortedfirst(wₛ, rand()*last(wₛ))
            μ += data[row] + randn()*sigma[row]
        end
        means[h] = μ/nd
    end
    means
end

## resample elemental compositions or ratios with weights and option of population abundances.
    #Replacement for simpleweightedresample
function weightedresample(df::DataFrame,vars::Vector{Symbol};
                            wt::Vector{Symbol},
                            wt_pop=nothing,
                            n::Int64=10_000,
                            μ::Bool=true,
                            ratio::Symbol=:none,
                            fraction::Bool=true)

    # Run weighted resampling through each category in 'vars' individually
    # to account for heterogeneous measurements
    d = copy(df)
    nᵣₛ = similar(vars,Int64)
    bsr_out = Array{Float64}(undef,n,length(vars))

# element for-loop
    @inbounds for k = 1:length(vars)
        j = vars[k]

#### Ratio resample option here
        if ratio == :none
            indᵣₛ = .!isnan.(d[:,j])
        else
            indᵣₛ = .!isnan.(d[:,j] .* d[:,ratio])
        end

        nᵣₛ[k] = sum(indᵣₛ) # number of non-NaN measurements
        meas = d[indᵣₛ,j]
        σₘ   = d[indᵣₛ,String(j)*"_err"]

        w    = Array{Float64}(undef,nᵣₛ[k],length(wt))
            # Calculate weights
        @inbounds for x = 1:length(wt)
            if wt_pop == nothing    # flat population profile
                w[:,x] = weights(d[indᵣₛ,wt[x]])
            elseif length(wt_pop[x]) > 0 # rough population profile
                w[:,x] = weights(d[indᵣₛ,wt[x]],wt_pop[x])
            else    # flat population profile but other rough profiles
                w[:,x] = weights(d[indᵣₛ,wt[x]])
            end
        end

        Πw =  vec(prod(w,dims=2))

        if ratio == :none
            μ ? bsr_out[:,k] = bootstrapmean(n,meas,σₘ,Πw) :
                bsr_out[:,k] = bootstrapresample(n,meas,σₘ,Πw)
        else
            bsr_out[:,k] =
                bootstrapratio(n,meas,d[indᵣₛ,ratio],Πw,
                            num_sigma = σₘ,
                            denom_sigma = d[indᵣₛ,String(ratio)*"_err"],
                            fraction=fraction,
                            means=μ)
        end
    end
    return (bsr_out, nᵣₛ)
end


## FRACTION AND RATIO RESAMPLING SCHEME
## UNDER CONSTRUCTION ??
# This assumes a gaussian analytical distribution

function bootstrapratio(nrs::Int64,num::AbstractVector,denom::AbstractVector,weights::AbstractVector=ones(Float64,size(num));
        num_sigma::AbstractVector=zeros(size(num)),
        denom_sigma::AbstractVector=zeros(size(num)),
        fraction::Bool=true,
        means::Bool = true
    )

    resampled = similar(num,nrs)
    wₛ = cumsum(weights)

    if means # Resample the full dataset `nrs` times with replacement & return means
        rsds = similar(num) # = `ReSampled DataSet`
        if fraction
            @inbounds for i = 1:nrs
                @inbounds for j=1:length(num)
                    row = searchsortedfirst(wₛ, rand()*last(wₛ))
                    x = num[row] + randn()*num_sigma[row]
                    rsds[j] =  x / (x+ denom[row] + randn() * denom_sigma[row])
                end
                μᵣₛ = mean(rsds)
                resampled[i] = μᵣₛ / (1-μᵣₛ)
            end
        else #calculate as ratios
            @inbounds for i = 1:nrs
                @inbounds for j=1:length(num)
                    row = searchsortedfirst(wₛ, rand()*last(wₛ))
                    rsds[j] =  (num[row] + randn()*num_sigma[row]) /
                        (denom[row] + randn() * denom_sigma[row])
                end
                resampled[i] = vmean(rsds)
            end
        end
    else # Resample `nrs` compositions from dataset.
    #rows = similar(num, Int64, nrs)
        if fraction
            @inbounds for i = 1:nrs
                row = searchsortedfirst(wₛ, rand()*last(wₛ))
                #rows[i] = row
                x = num[row] + randn()*num_sigma[row]
                resampled[i] =  x / (x + denom[row] + randn() * denom_sigma[row])
            end
        else
            @inbounds for i = 1:nrs
                row = searchsortedfirst(wₛ, rand()*last(wₛ))
                #rows[i] = row
                resampled[i] =  (num[row] + randn()*num_sigma[row]) /
                    (denom[row] + randn() * denom_sigma[row])
            end
        end
    end
    return resampled
end

## Resample from a dataframe for multiple types, with the "type" column given as type=:column
    # e.g. resample for each petrologic type with `type=:pet_type`
function wrstype(df::DataFrame,vars::Vector{Symbol};
                            type::Symbol,
                            wt::Vector{Symbol},
                            wt_pop=nothing,
                            n::Int64=10_000,
                            μ::Bool=true,
                            ratio::Symbol=:none,
                            fraction::Bool=true)
    d=copy(df)

    types = sort(unique(d[:,type]))

    rs_out = Array{Float64}(undef,n,length(vars),length(types))
    nₜ = Array{Int64}(undef,length(vars),length(types))

    for i = 1:length(types)
        d_wrs = d[d[:,type] .== types[i],:]
        ( rs_out[:,:,i], nₜ[:,i] ) = weightedresample(d_wrs, vars,wt=wt,wt_pop=wt_pop,n=n,μ=μ,ratio=ratio,fraction=fraction)
    end
    return (types,rs_out,nₜ)
end

#= Weighted Resampling Function

function simpleweightedresample(d::DataFrame,vars::Vector;n::Int64,μ::Bool=false)

# [?] Boolean to output metadata
# [x] Multiple weighting schemes
    #~ Fixed with "weights", with which multiple weights can be calculated and combined

# Run weighted resampling through each category in 'vars' individually
# to account for heterogeneous measurements
    nⱼ = similar(vars,Int64)
    bsr_out = Array{Float64}(undef,n,length(vars))
    for k = 1:length(vars)
        j = vars[k]
        indⱼ = .!isnan.(d[:,j])
        nⱼ[k] = sum(indⱼ) # number of non-NaN measurements

#################
        name = d[indⱼ,:NAME] # Hardwired to weight only by :NAME column.
#################
        meas = d[indⱼ,j]
        σₘ   = d[indⱼ,String(j)*"_err"]
        w    = Vector{Float64}(undef,nⱼ[k])

# Step 2 identify and count discreet samples (e.g. specific meteorites)
        samples = unique(name)

        for i in samples
            indᵢ = (name .== i)
            nᵢ = count(indᵢ)    # Calculate counts of iᵗʰ sample
            w[indᵢ] .= 1 / nᵢ # Calculate weight from counts
        end

        μ ? bsr_out[:,k] = bootstrapmean(n,meas,σₘ,w) :
            bsr_out[:,k] = bootstrapresample(n,meas,σₘ,w)

    end
    return (bsr_out, nⱼ)
end
=#

