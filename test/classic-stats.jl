wm = SolarChem.weightedmean([3.4, 4.2, 5.1], [.4, .6, .7])
@test wm[1] ≈ 3.914212548015365
@test wm[2] ≈ 0.3005756321250956

uwm = SolarChem.unweightedmean([3.4, 4.2, 5.1], [.4, .6, .7])
@test uwm[1] == mean([3.4, 4.2, 5.1])
@test uwm[2] ≈ 1.004987562112089