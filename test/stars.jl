## solartwins

sttest = solartwins()
nogce = solartwins(gce=false)
@test sttest.Dy[79] ≈ 0.060058478
@test sttest.sDy[79] ≈ 0.019478009
@test sttest.Sc[1] ≈ (-0.011278679+ -0.011278679)/2
@test sttest.sSc[1] ≈ sqrt(0.008092654^2 + 0.013057222^2)

@test nogce.sDy[79] == 0.019
@test nogce.Dy[79] == 0.074
@test nogce.Sc[1] ≈ -0.019
@test nogce.sSc[1] ≈ sqrt(0.013^2 + 0.008^2)