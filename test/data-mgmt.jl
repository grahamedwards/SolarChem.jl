## Testing data-mgmt.jl functions
    # estimateuncertainty
    # calcweights
    
Na = [8.7, NaN, 8.0, 7.4, 8.2, NaN]
Ca = [17.6, 13.9, 14.9, NaN, 15.1, 15.6]
Mg = [36.1, 35.9, 35.8, NaN, 36.4, 36.1]
Fe = 18ones(40)
sFe = vcat(fill(.08*18,20), fill(1.8,17),[NaN,NaN,35])

d = (; Na, sNa = .6Na, Ca, sCa = fill(NaN,length(Ca)), Mg, Fe, sFe)


sigpct = 4
dd = estimateuncertainty(d, sigpct)
@test dd.sMg[1] == .01sigpct*dd.Mg[1] # create sMg 
@test dd.sNa[1] == .01sigpct*dd.Na[1] # correct uncs > maxpctunc
isnan(last(dd.sNa)) # NaN unc for NaN measurements.
# calculate mean unc and replacement.
    @test .08 < last(dd.sFe) / last(dd.Fe) < .09
    @test !isnan(dd.sFe[end-1])
    
v = ["a", "b", "a", "c",  "b", "a"]

@test calcweights(v) ≈ [1/3, 0.5, 0.1/3, 1.0, 0.5, 1/3]

@test calcweights(v, weights = Dict("a" => 1/3, "b" => 1/2, "c" => 1.0)) == ones(6)






