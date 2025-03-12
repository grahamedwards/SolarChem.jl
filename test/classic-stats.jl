wm = SolarChem.weightedmean([3.4, 4.2, 5.1], [.4, .6, .7])
@test wm[1] ≈ 3.914212548015365
@test wm[2] ≈ 0.3005756321250956

uwm_std = SolarChem.unweightedmean([3.4, 4.2, 5.1], [.4, .6, .7], sem=false)
uwm_sem = SolarChem.unweightedmean([3.4, 4.2, 5.1], [.4, .6, .7])
@test uwm_std[1] == uwm_sem[1] == mean([3.4, 4.2, 5.1])
@test uwm_std[2] ≈ 1.004987562112089
@test uwm_sem[2] ≈ uwm_std[2]/sqrt(3)