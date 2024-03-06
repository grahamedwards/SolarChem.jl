

x, xsig, xwts  = [1.8, 2.1, 2.6, 2.3], [.3, .2, .2, .3], [0.5, 1, 0.5, 1]

@test bsresample(3,x,xsig,w=xwts,rng=StableRNG(4567)) == [2.588272554749283, 1.5156891982328184, 2.1409336364526705]


@test bsmean(4,x,xsig,w=xwts, rng=StableRNG(4567)) == [2.2027592881961953, 2.3916504384860113, 2.195010862905976, 2.2217297502067046]

d = (name = ["b", "b", "a", "c", "b", "a"],  Na = [8.7, NaN, 8.0, 7.4, 8.2, NaN], Ca = [17.6, 13.9, 14.9, NaN, 15.1, 15.6], Mg = [36.1, 35.9, 35.8, NaN, 36.4, 36.1], sNa = [0.348, NaN, 0.32, 0.296, 0.328, NaN], sCa = [0.704, 0.556, 0.596, NaN, 0.604, 0.624], sMg = [1.444, 1.436, 1.432, NaN, 1.456, 1.444])

x = bootstrapelements(2,d,(:Na,:Mg), rng=StableRNG(1))
    @test x.Na[2] ≈ 8.293674810925634
    @test x.Mg[2] ≈  35.94359320366079

x = bootstrapelements(2,d,(:Na,:Mg), resamplemeans=false, rng=StableRNG(1))
    @test x.Na[2] ≈ 7.150312928453221
    @test x.Mg[2] ≈ 37.38362424031754


x = bootstrapelements(2,d,(:Na,:Mg), weighted=false, rng=StableRNG(1)) 
    @test x.Na[2] ≈ 8.490977284313256
    @test x.Mg[2] ≈ 35.962586786415926


x = bootstrapratios(2,d,(:Na,:Mg), :Ca, rng=StableRNG(1)) 
    @test x.Na[2] ≈ 0.6613934042322102
    @test x.Mg[2] ≈ 1.8953207799735892


x = bootstrapratios(2,d,(:Na,:Mg), :Ca, resamplemeans=false, rng=StableRNG(1)) 
    @test x.Na[2] ≈ 0.26698078243809087
    @test x.Mg[2] ≈ 3.394220266807608

x = bootstrapratios(2,d,(:Na,:Mg), :Ca, weighted=false, rng=StableRNG(1)) 
    @test x.Na[2] ≈ 0.6652394491425418
    @test x.Mg[2] ≈ 1.94781996814001

