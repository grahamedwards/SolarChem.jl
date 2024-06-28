"""

    solartwins(; gce=true)

Returns a named tuple containing stellar parameters and solar-normalized compositions ([X/Fe], in dex) either corected for galactic chemical evolution (`gce=true`) or left uncorrected (`gce=false`). All corresponding uncertainty fields are preceded with `s`, e.g. `:Ca` and `:sCa` or `logg` and `slogg`.

Data is compiled from Bedell + 2018 (*ApJ*, doi:[10.3847/1538-4357/aad908](https://doi.org/10.3847/1538-4357/aad908)). Where multiple species are reported (e.g. CI-CH, ScI-ScII, TiI-TiII, CrI-CrII), this reports the mean and quadrature-propagated uncertainties of the two species.

"""
function solartwins(;gce::Bool=true)

    x = DelimitedFiles.readdlm(string(@__DIR__,"/../data/Bedell2018-solar-twins.csv"),',')

    view(x,12,21:24) .= NaN

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