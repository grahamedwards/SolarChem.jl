## Testing data-mgmt.jl functions
    # estimateuncertainty
    # calcweights

Na = [8.7, NaN, 8.0, 7.4, 8.2, NaN]
Ca = [17.6, 13.9, 14.9, NaN, 15.1, 15.6]
Mg = [36.1, 35.9, 35.8, NaN, 36.4, 36.1]

d = (; Na, sNa = .6Na, Ca, sCa = fill(NaN,length(Ca)), Mg)

dd = estimateuncertainty(d, 4)
dd.sMg[1] == 0.04dd.Mg[1]
dd.sNa[1] == 0.04dd.Na[1]
isnan(last(dd.sNa))



    
v = ["a", "b", "a", "c",  "b", "a"]

@test calcweights(v) ≈ [0.3333333333333333, 0.5, 0.3333333333333333, 1.0, 0.5, 0.3333333333333333]

@test calcweights(v, weights = Dict("a" => 1/3, "b" => 1/2, "c" => 1.0)) == ones(6)






