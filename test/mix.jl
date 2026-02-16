@test SolarChem.normpdf(2.2, 2., 0.4) ≈ 0.8801633169107486

@test exp(SolarChem.fastll(2.2, 2., 0.4))/(sqrt(2π)*.4) ≈ SolarChem.normpdf(2.2, 2., 0.4)

# lldist interior

@test -1.5 > SolarChem.lldist(2.,1., 3.,1.,n=100) > -1.6 

@test abs(SolarChem.lldist(2.,1., 3.,1.,n=1000) - SolarChem.lldist(2.,1., 3.,1.)) < 0.0001 

# lldist outer

tm,ts = [2., 2., 2.], [1., 1., 1.]

@test SolarChem.lldist(tm, ts, 3., 1.,n=100) ≈ 3*SolarChem.lldist(2.,1., 3.,1.,n=100)

@test SolarChem.lldist(tm, ts, 3., NaN) == SolarChem.lldist(tm, ts, NaN, 1.) == 0


pt = jt = Fractions(0.8, 0.8)


@test SolarChem.jump(pt, jt, rng=StableRNG(11))[2] == :outer
@test SolarChem.jump(pt, jt, rng=StableRNG(11))[1].outer  ==  SolarChem.jump(pt, jt, rng=StableRNG(11))[3] + pt.outer


@test SolarChem.solarlogmix(Composition(.003,1e-4), Composition(.04, 1e-3), Composition(0.3, .03), pt)[1] ≈ log10((pt.sun * (0.003 * (1-pt.outer) + .04 * pt.outer) + 0.3)/(pt.sun+1)/0.3)[1]

@test SolarChem.solarlogmix(Composition(.003,1e-4), Composition(.04, 1e-3), Composition(0.3, .03), pt)[2] > SolarChem.solarlogmix(Composition(.003,1e-4), Composition(.04, 1e-3), Composition(0.3, .03), pt,solarunc=false)[2]



# solarmixmetopolis test later...

