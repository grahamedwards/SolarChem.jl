
"""

    weightedmean(m, s)

Given values `m` (e.g. measurement means) with absolute 1σ uncertainties `s`, returns the variance-weighted mean and 1σ uncertainty. 

``\\mu = \\frac{\\sum^n_{i=1} x_i~\\sigma_i^{-2}}{\\sum^n_{i=1}\\sigma_i^{-2}}``

``\\sigma = \\sqrt{\\frac{1}{\\sum^n_{i=1}\\sigma_i^{-2}}}``

"""
function weightedmean(m::Vector, s::Vector) 
    @assert length(m) == length(s)
    n = length(m)
    Sw = Sxw = 0.0
    @inbounds @simd for i = eachindex(m)
        si = s[i]
        wi = 1.0 / (si*si)
        Sw += wi
        Sxw += m[i] * wi 
    end

    ( Sxw/Sw , sqrt(1/Sw))  
end

"""

    unweightedmean(m, s; sem=false)

Given values `m` (e.g. measurement means) with absolute 1σ uncertainty `s`, returns the unweighted mean and 1σ uncertainty estimated by addition in quadrature. 

``\\sigma = \\sqrt{\\sum^n_{i=1} s_i^2}``

If `sem`=`true`, the 1σ uncertainty is the standard error of the mean (``\\sigma_{\\bar x} = \\sigma/\\sqrt{n}``).

"""
function unweightedmean(m::Vector, s::Vector; sem::Bool=true)
    @assert length(m) === length(s)
    Statistics.mean(m), sqrt(sum(s.^2)/ifelse(sem,length(s),1.0))
end