## Testing data-mgmt.jl functions
    # estimateuncertainty
    # calcweights
    # trimnans
    # pulltopic
    # pulltype
    # pullgroup
    # exclude
    
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

sFe .= NaN
@test estimateuncertainty(d, sigpct).sFe[end] == .04last(Fe)
    
v = ["a", "b", "a", "c",  "b", "a"]

@test calcweights(v) ≈ [1/3, 0.5, 1/3, 1.0, 0.5, 1/3]
@test calcweights(v, weights = Dict("a" => 1/3, "b" => 1/2, "c" => 1.0)) == ones(6)

name = comment = citation = dataset = v

group = [(:H,), (:H, :L), (:H,), (:CV,), (:H,), (:L,)] 
type = [(3,), (3,), (5,6), (3,), (7,), (5,)]

d = estimateuncertainty((; name, type, group, comment, citation, dataset, Na, Ca, Mg), sigpct)

@test trimnans(d,:Na).name == ["a", "a", "c", "b"]
@test trimnans(d,:Na).Na == Na[.!isnan.(Na)]
@test trimnans(d,(:Na,:Mg)).Mg == Mg[.!isnan.(Na) .* .!isnan.(Mg)]


#@test @suppress prod(isnan.(trimextremes(d).Ca))
#@test @suppress count(isnan.(trimextremes(d, min=14, max = 20).Ca)) == 2

@test pulltopic(d, :comment, "c").comment == ["c"]
@test length(pulltopic(d, :comment, "c").group) == length(pulltopic(d, :comment, "c").Na)
@test pulltopic(d, :comment, ("b","c")).comment == ["b", "c", "b"]
d.comment[1] = "bc"
@test pulltopic(d, :comment, "c").comment == ["bc", "c"]
@test pulltopic(d, :comment, "c", exactmatch=true).comment == ["c"]


@test length(exclude(d,:comment,"a").name) == length(exclude(d,:comment,"a").Na)
@test exclude(d,:comment,"a").comment == ["bc", "b", "c", "b"]
@test exclude(d,:comment,("a","b")).comment == ["c"]


dheat =deepcopy(d)
dheat.comment .= ["bc °C", "370 degrees C", "350 C", "355 C", "b", "a"]
@test excludeheated(dheat).name == ["b", "a"]


@test length(pullgroup(d,:H).group) == length(pullgroup(d,:H).Na)
@test pullgroup(d,:H).group == [(:H,), (:H, :L), (:H,), (:H,)] 
@test pullgroup(d,(:H,:L)).group == [(:H,), (:H, :L), (:H,), (:H,), (:L,)] 
@test pullgroup(d,(:H,:L), exactmatch=true).group == [(:H, :L)] 

@test length(pulltype(d,3).group) == length(pulltype(d,3).Na)
@test pulltype(d,3).type == fill((3,),3)
@test pulltype(d,(3,5)).type == [(3,), (3,), (5,6), (3,), (5,)] 
@test pulltype(d,(5,6), exactmatch=true).type == [(5,6)] 

@test countratios(d,:Na,:Ca) == 3
@test countratios(d,(:Na,:Mg),:Ca) == (Na = 3, Mg = 5, divisor = :Ca)
@test countmeasurements(d,:Na) == 4
@test countmeasurements(d,(:Na,:Mg)) == (Na = 4, Mg = 5)