Requires.@require CairoMakie="13f3f980-e62b-5c42-98c6-ff1f3baf88f0" @eval using .CairoMakie
Requires.@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a"  @eval using .WGLMakie
Requires.@require WGLMakie="276b4fcb-3e11-5398-bf8b-a0c2d153d008" @eval using .WGLMakie

export cleanhist, onehist, histpanels

"""
```julia
interleave(a)
```
Interleave the sequential values of a vector `a` so that each value occurs twice in a vector of `length(a)*2`

### Example

    julia> interleave([1,2,3])
    6-element Vector{Int64}:
    1
    1
    2
    2
    3
    3

"""
function interleave(a::AbstractVector{T}) where T
    c = Vector{T}(undef, 2*(length(a)))
    interleave!(c,a)
    return c
end

"""
```julia
interleave!(c,a)
```
In-place version of [`interleave`](@ref), that fills `c` with pairs of each index in `a`.
"""
function interleave!(c::Vector{T}, a::AbstractVector{T}) where T
    @assert length(c) == 2*length(a)
    i = 0
    @inbounds for x in 1:length(a)
        c[i += 1] = a[x]
        c[i += 1] = a[x]
    end
    return c
end


"""

```julia
binweave(a)
```

Interleave the sequential values of `a` (`<:AbstractVector`) using the first and last values only once. To be used with [`interleave`](@ref) to make plotting-ready histograms.

### Example

    julia> binweave([1,2,3])
    4-element Vector{Int64}:
    1
    2
    2
    3


"""
function binweave(a::AbstractVector{T}) where T
    c = Vector{T}(undef, 2*(length(a)-1))
    binweave!(c,a)
    return c
end


"""
```julia
binweave!(c,a)
```
In-place version of [`binweave`](@ref), that overwrites `c` with woven bin edges.

"""
function binweave!(c::Vector{T}, a::AbstractVector{T}) where T
    @assert length(c) == 2*(length(a)-1)
    i = 0
    @inbounds for x in 1:length(a)-1
        c[i += 1] = a[x]
        c[i += 1] = a[x+1]
    end
    c 
end



function quickhist(x::AbstractArray, binedges::AbstractRange)
    h = zeros(Int, length(binedges)-1)
    binmin, binstep, binmax = first(binedges), step(binedges), last(binedges)
    @inbounds for i in eachindex(x)
        xi = x[i]
        binmin < xi <= binmax  || continue
        h[1+floor(Int, (xi - binmin)/binstep)] += 1
    end
    h
end

"""

```julia
cleanhist(x; nbins=32, scooch_nbins=2)
```
Calculates a histogram with extra (0 count) bins to buffer the edges and make it look nice and clean. 🧼

Optionally specify the number of histogram `bins` (default: `2⁵` bins) and the number of buffering bins `scooch`. (Total bins = `nbins + scoochbins`)

Returns a `NamedTuple` with `x` and `y` values of histogram.

"""
function cleanhist(x::Vector{<:Number}; bins::Int=32, scooch::Int=2)
    xmin,xmax = extrema(x)
    x_scooch = scooch*(xmax-xmin)/(bins)
    binedges = LinRange(xmin-x_scooch, xmax+x_scooch, bins+2*scooch+1)
    y = quickhist(x, binedges) ./ (length(x)*step(binedges))
    return (x=SolarChem.binweave(binedges), y=SolarChem.interleave(y))
end


function histpanels(data_in::NamedTuple; els::Tuple{Vararg{Symbol}}=(), cols::Int=0, bins::Int=32, labelsuffix=" (m/m)", figsize=(800,600), darkmode::Bool=false)

    els = ifelse(isempty(els), keys(data_in), els)
    nels = length(els)
    cols = ifelse(iszero(cols), ceil(Int,sqrt(nels)), cols)
    rows = ceil(Int,nels/cols)

    f= Figure(size=(figsize), backgroundcolor=ifelse(darkmode,:transparent,:white))
    for j in CartesianIndices((rows,cols))
        i = j[1] + rows * (j[2]-1) #calculate index in els
        if i <= nels
            el = els[i] 
            onehist(data_in[el], el, f=f[j[1],j[2]], bins=bins, labelsuffix=labelsuffix, darkmode=darkmode)
        end
    end
    f
end

function onehist(data_in::Vector, el::Symbol; f= Figure(), bins::Int=32, labelsuffix=" (g/g)", darkmode::Bool=false)
    x=data_in
# Convert
    pltclr = fillcolor = ifelse(darkmode,:white,:black)

    h = SolarChem.cleanhist(x, bins=bins,scooch=2)

    ax = Axis(f[1,1], xlabel=string(el,labelsuffix),bottomspinecolor=pltclr,xtickcolor=pltclr,xticklabelcolor=pltclr, xlabelcolor=pltclr,backgroundcolor=ifelse(darkmode,:transparent,:white),
    xgridvisible=false,ygridvisible=false,yticklabelsvisible=false,yticksvisible=false,rightspinevisible=false,leftspinevisible=false,topspinevisible=false,)

    band!(ax,h.x,h.y,zero(h.y), color=(fillcolor,0.1))
    lines!(ax,h.x,h.y, color=pltclr,linewdith=2,)
    f
end

