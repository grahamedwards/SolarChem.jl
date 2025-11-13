testmeas = (; Ca = [2., 2.8, 1.8, 2.3, 1.7], Fe = [8.6, 7.2, 7.8, 7.6, 8.6], Mg = [3.3, 2.5, 2.4, 3.5, 2.4] )

testratio = (; Ca = [ .25, .34, .21, .28], Mg=[.4, .37, .28, .42, .25], divisor=:Fe)

ntest = (; Ca= 6, Mg = 3, divisor=:Fe)

@test ratiosummary(testmeas,testratio)[1] == collect(keys(testratio)[1:end-1]) # test the returned numerator identifiers.

@test ratiosummary(testmeas,testratio)[2] ≈ [mean(testratio.Ca) std(testratio.Ca) ; mean(testratio.Mg) std(testratio.Mg)] # test ratio mean/std calc.

@test ratiosummary(testmeas,testratio)[2] ≈ ratiosummary(testmeas,testratio,ratiocounts=ntest, minratios=2)[2] # All ratioscounts >  minratios

@test ratiosummary(testmeas,testratio,ratiocounts=ntest, minratios=4)[2][2,1] == mean(testmeas.Mg)/mean(testmeas.Fe) # Test mean calculation from measurements. 

@test ratiosummary(testmeas,testratio,ratiocounts=ntest, minratios=4)[2][2,2] == mean(testmeas.Mg)/mean(testmeas.Fe) * sqrt((std(testmeas.Mg)/mean(testmeas.Mg))^2 + (std(testmeas.Fe)/mean(testmeas.Fe))^2) # Test std calculation from measurements. 

# test ratio calculation of all NaN fields. 
trnan = deepcopy(testratio)
trnan.Mg .= NaN

@test ratiosummary(testmeas,testratio,ratiocounts=ntest, minratios=4)[2] ≈ ratiosummary(testmeas,trnan)[2] 

# test NamedTuple output

@test ratiosummary(testmeas,testratio, NTout=true) isa NamedTuple