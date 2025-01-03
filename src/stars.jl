"""

    solartwins(; gce=true, includesun=false)

Returns a named tuple containing stellar parameters and solar-normalized compositions ([X/Fe], in dex) either corected for galactic chemical evolution (`gce=true`) or left uncorrected (`gce=false`). All corresponding uncertainty fields are preceded with `s`, e.g. `:Ca` and `:sCa` or `logg` and `slogg`. By default excludes the Sun (improves mixing model efficiency), but you can include solar data (mostly 0 dex by definition) with includesun=`true`.

Data is compiled from Bedell + 2018 (*ApJ*, doi:[10.3847/1538-4357/aad908](https://doi.org/10.3847/1538-4357/aad908)). Where multiple species are reported (e.g. CI-CH, ScI-ScII, TiI-TiII, CrI-CrII), this reports the mean and quadrature-propagated uncertainties of the two species.

"""
function solartwins(;gce::Bool=true, includesun::Bool=false)

    x = DelimitedFiles.readdlm(string(@__DIR__,"/../data/Bedell2018-solar-twins.csv"),',')

    x = includesun ? x : x[1:end-1, :]

    @inbounds @simd ivdep for i=eachindex(x)
        xi = x[i]
        x[i] = ifelse(isempty(xi),NaN,xi)
    end 

    heads = view(x,1,:)

    elsym = allsolar()

    name = Symbol.(x[2:end,1])
    thickdisk = Bool[ifelse(xx =="True",true,false) for xx in x[2:end,12]]
    starnumbers = NamedTuple{(:teff, :steff, :logg, :slogg, :feh, :sfeh, :age, :sage, :mass, :smass)}(float.(x[2:end,i]) for i=2:11)

    starstuff=(;name, starnumbers..., thickdisk)

    y = NamedTuple{(Symbol.(heads[13:end])...,)}(float.(x[2:end,i]) for i= 13:size(x,2));

    n,v,vgce = (),(),()

    @inbounds for i = eachindex(elsym)
        el = elsym[i]
        if el == :Fe
        elseif el == :C 
            C = @. (y.CI + y.CH) / 2
            sC = @. sqrt(y.CI_err^2 + y.CH_err^2)
            Cg =  @. (y.CI_gce + y.CH_gce) / 2
            sCg = @. sqrt(y.CI_gce_err^2 + y.CH_gce_err^2)
            n,v,vgce = (n...,:C,:sC),(v...,C,sC),(vgce...,Cg,sCg)

        elseif el ∈ (:Sc, :Ti, :Cr)
            elv = @. (y[Symbol(el,:I)] + y[Symbol(el,:II)]) / 2
            selv= @. sqrt(y[Symbol(el,:I_err)]^2 + y[Symbol(el,:II_err)]^2)
            gelv = @. (y[Symbol(el,:I_gce)] + y[Symbol(el,:II_gce)]) / 2
            gselv= @. sqrt(y[Symbol(el,:I_gce_err)]^2 + y[Symbol(el,:II_gce_err)]^2)
            n,v,vgce = (n...,el,Symbol(:s,el)),(v...,elv,selv),(vgce...,gelv,gselv)
        else
            spI, spII = Symbol.(el,(:I,:II))
            sp = ifelse(spI ∈ keys(y), spI, spII)
            n = (n...,el,Symbol(:s,el))
            v = (v...,y[sp], y[Symbol(sp,:_err)])
            vgce = (vgce...,y[Symbol(sp,:_gce)], y[Symbol(sp,:_gce_err)])
        end
    end

    comps = gce ? NamedTuple{n}(vgce) : NamedTuple{n}(v)

    (; starstuff..., comps...)
end


lingce(m::Number,b::Number,age::Number, x::Number) = x - (m*age+b)

function quadgce(aa::Number, bb::Number, age::Number, x::Number)
    da = age-4.6
    x - da*(aa + bb * da)
end


"""

    sun(divisor=nothing; logspace=true)

Returns a NamedTuple containing the solar photosphere composition, as reported by Asplund+ 2021 [*A&A*, doi:[10.1051/0004-6361/202140445](https://doi.org/10.1051/0004-6361/202140445)]. Optionally provide an element name for `divisor`as a Symbol, to return the `divisor`-normalized ratios with uncertainties propagated in quadrature. `logspace=false` converts values and uncetainties to linear-space ratios. 

"""
function sun(divisor::Union{Symbol,Nothing}=nothing; logspace::Bool=true)

    elements = (:H, :He, :Li, :Be, :B, :C, :N, :O, :F, :Ne, :Na, :Mg, :Al, :Si, :P, :S, :Cl, :Ar, :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Kr, :Rb, :Sr, :Y, :Zr, :Nb, :Mo, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn, :Sb, :Te, :I, :Xe, :Cs, :Ba, :La, :Ce, :Pr, :Nd, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, :Hf, :Ta, :W, :Re, :Os, :Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Th, :U)

    isnothing(divisor) || @assert divisor ∈ elements

    x = ([12.0, 0.0], [10.914, 0.013], [0.96, 0.06], [1.38, 0.09], [2.7, 0.2], [8.46, 0.04], [7.83, 0.07], [8.69, 0.04], [4.4, 0.25], [8.06, 0.05], [6.22, 0.03], [7.55, 0.03], [6.43, 0.03], [7.51, 0.03], [5.41, 0.03], [7.12, 0.03], [5.31, 0.2], [6.38, 0.1], [5.07, 0.03], [6.3, 0.03], [3.14, 0.04], [4.97, 0.05], [3.9, 0.08], [5.62, 0.04], [5.42, 0.06], [7.46, 0.04], [4.94, 0.05], [6.2, 0.04], [4.18, 0.05], [4.56, 0.05], [3.02, 0.05], [3.62, 0.1], [NaN, NaN], [NaN, NaN], [NaN, NaN], [3.12, 0.1], [2.32, 0.08], [2.83, 0.06], [2.21, 0.05], [2.59, 0.04], [1.47, 0.06], [1.88, 0.09], [1.75, 0.08], [0.78, 0.11], [1.57, 0.1], [0.96, 0.1], [NaN, NaN], [0.8, 0.2], [2.02, 0.1], [NaN, NaN], [NaN, NaN], [NaN, NaN], [2.22, 0.05], [NaN, NaN], [2.27, 0.05], [1.11, 0.04], [1.58, 0.04], [0.75, 0.05], [1.42, 0.04], [0.95, 0.04], [0.52, 0.04], [1.08, 0.04], [0.31, 0.1], [1.1, 0.04], [0.48, 0.11], [0.93, 0.05], [0.11, 0.04], [0.85, 0.11], [0.1, 0.09], [0.85, 0.05], [NaN, NaN], [0.79, 0.11], [NaN, NaN], [1.35, 0.12], [NaN, NaN], [NaN, NaN], [0.91, 0.12], [NaN, NaN], [0.92, 0.17], [1.95, 0.08], [NaN, NaN], [0.03, 0.1], [NaN, NaN])

    xx= ()

    if isnothing(divisor)
        @inbounds for i = eachindex(x)
            xi = x[i]
            xi[1] -= 12.0
            xi1 = 10^xi[1]
            xi2 = log(10) * xi1 * xi[2]
            xx= (xx..., Composition(xi1, xi2))
        end
    else
        id = findfirst( elements .== divisor)
        d, ud2 = x[id][1], x[id][2]^2
    
        @inbounds for i = eachindex(x)
            xi = x[i]
            xm  = xi[1] - d
            xs = sqrt(xi[2]*xi[2] + ud2)
            xms =  logspace ? (xm,xs) : (10^xm, log(10) * xs * 10^xm )
            xx = (xx..., Composition(xms...))
        end 
    end
    (; zip(elements, xx)...)
end

