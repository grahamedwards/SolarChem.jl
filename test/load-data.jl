## Groups 
@test SolarChem.assigngroup("EH-IMP MELT") == (:EH,)
@test SolarChem.assigngroup("EH4") == (:EH,)
@test SolarChem.assigngroup("EL6/7") == (:EL,)
@test SolarChem.assigngroup("E6") == (:E,)
@test SolarChem.assigngroup("ENSTATITE METEORITE UNGR") == (:E,)
@test SolarChem.assigngroup("ENST ACHON") == ()

@test SolarChem.assigngroup("H5") == (:H,)
@test SolarChem.assigngroup("H-IMP MELT") == (:H,)
@test SolarChem.assigngroup("L5") == (:L,)
@test SolarChem.assigngroup("LL6") == (:LL,)
@test SolarChem.assigngroup("H/L3.9") == (:H,:L)
@test SolarChem.assigngroup("L/LL4") == (:L,:LL)
@test SolarChem.assigngroup("L(LL)3.05 CHONDRITE") == (:L,:LL)

@test SolarChem.assigngroup("CI1") == (:CI,)
@test SolarChem.assigngroup("CV2") == (:CV,)
@test SolarChem.assigngroup("CM1/2") == (:CM,)
@test SolarChem.assigngroup("EUCRITE-CM") == ()
@test SolarChem.assigngroup("CR2") == (:CR,)
@test SolarChem.assigngroup("CH3") == (:CH,)
@test SolarChem.assigngroup("CH3") == (:CH,)
@test SolarChem.assigngroup("CBB") == (:CB,)
@test SolarChem.assigngroup("CBA") == (:CB,)
@test SolarChem.assigngroup("CK4/5") == (:CK,)
@test SolarChem.assigngroup("CL4") == (:CL,)
@test SolarChem.assigngroup("CO3.0") == (:CO,)
@test SolarChem.assigngroup("C1/2-UNG") == (:Cung,)
@test SolarChem.assigngroup("LUNAR-FELDSPATHIC BRECCIA") == ()


## Types

@test SolarChem.assigntype("H5") == (5,)
@test SolarChem.assigntype("LL3.15") == (3,)
@test SolarChem.assigntype("L/LL3.5") == (3,)
@test SolarChem.assigntype("LL5/6") == (5, 6)
@test SolarChem.assigntype("L3.3-3.6 CHONDRITE") == (3,)
@test SolarChem.assigntype("R3.8-5") == (3, 4, 5)

@test SolarChem.assigntype("H~6") == (6,)
@test SolarChem.assigntype("H-IMP") == (7,)
@test SolarChem.assigntype("CM-AN") == (0,)
@test SolarChem.assigntype("H-AN") == (0,)
@test SolarChem.assigntype("H3-AN") == (3,)


@test SolarChem.unitconversionfactor.("Ca" .* ("wt%", "%", "%m/m", "mg/g", "μg/g", "ppm", "ppb", "ng/g", "pg/g")) == (0.01, 0.01, 0.01, 0.001, 1.0e-6, 1.0e-6, 1.0e-9, 1.0e-9, 1.0e-12)
@test @silence isnan(SolarChem.unitconversionfactor("boop"))

@test SolarChem.readvalue(1) == 1.
@test SolarChem.readvalue(2.3) == 2.3
@test isnan(SolarChem.readvalue(""))
@test SolarChem.readvalue("1.2e8") == 1.2e8

@test SolarChem.unitconversionfactor("Ca (wt%)") === 0.01
@test SolarChem.unitconversionfactor("Ca (mg/g)") === 0.001
@test SolarChem.unitconversionfactor("Ca (ppm)") === 1e-6
@test SolarChem.unitconversionfactor("Ca (pg/g)") === 1e-12


#= Things to test in loadastromatdata()
    1. Compile headers: sample, taxonName, analysisComment, citation_url, dataset_url, a duplicate dataset_url, and a throwaway ignoreMe
    2. Oxides: CaO, FeO, and Fe2O3 columns.
    3. lab, method, and comment columns. 
    4. Appropriate scaling of a few units.
    5. One relative uncertainty, one S absolute, and one 2S uncertainty.
    6. One no uncertainty column. 
=#
testromat = loadastromatdata(file="testromat.csv")

@test keys(testromat) === (:name, :group, :type, :comment, :citation, :dataset, :Ca, :Fe, :sFe, :Cu, :sCu, :U, :sU)
@test testromat.name == ["Rock", "Em", "Sock", "Em", "Robots"]
@test unique(testromat.citation) == ["https://github.com/grahamedwards/SolarChem.jl"]
@test unique(testromat.dataset) == ["data-url"]
@test testromat.type == [(3, 4, 5, 6, 7), (5,), (4,), (5,), (3,)]
@test testromat.group == [(:LL,), (:H,), (:L,), (:H,), (:CV,)]
@test testromat.comment == ["dark matrix", "fusion crust", "on buttered toast", "ICP-MS", "science"]

@test 15.1 * 0.007146 < testromat.Ca[1] < 0.007148 * 15.1
@test isnan(testromat.Ca[3])
@test 13 * 0.007146 < testromat.Ca[5] < 0.007148 * 13

@test isapprox(0.01(22 * 0.77731 + 4.8 * 0.69943), testromat.Fe[1],atol=0.0001)
@test isapprox( 0.02 * 0.69943, testromat.Fe[4],atol=0.0001)
@test isapprox( 0.189 * 0.77731, testromat.Fe[3],atol=0.0001)
@test isapprox( 0.012 * 0.77731, testromat.sFe[3],atol=0.0001)
@test isapprox( sqrt((0.023 * 0.77731)^2 + (0.004 * 0.69943)^2), testromat.sFe[1],atol=0.0001)
@test isapprox( .1(.24 * 0.77731 +  0.05 * 0.69943), testromat.sFe[5],atol=0.0001)

@test testromat.Cu[[1,3,5]] ≈ [0.00395, 0.0743, 0.06]
@test testromat.sCu[[1,3]] ≈ [ 0.034365, 0.15603]
@test testromat.U[[1,4]] ≈ [1e-10, 2e-10]
@test testromat.sU[[1,4]] ≈ [2e-12,2.5e-11]

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