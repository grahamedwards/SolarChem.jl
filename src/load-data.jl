




"""

    SolarChem.periodictable

Tuple containing all elements (as `Symbol`s) in order of increasing Z, up to Z=96. A constant within the scope of the SolarChem module. 

"""
const periodictable = (
    :H, :He, # row 1
    :Li, :Be, :B, :C, :N, :O, :F, :Ne, # row 2
    :Na, :Mg, :Al, :Si, :P, :S, :Cl, :Ar, # row 3
    :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Kr, # row 4
    :Rb, :Sr, :Y, :Zr, :Nb, :Mo, :Tc, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn, :Sb, :Te, :I, :Xe, # row 5
    :Cs, :Ba, # row 6, groups 1-2
    :La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, # REEs (row 6, group 3)
    :Hf, :Ta, :W, :Re, :Os, :Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Po, :At, :Rn, # row 6, groups 4-18
    :Fr, :Ra, :Ac, :Th, :Pa, :U, :Np, :Pu, :Am, :Cm) #row 7, some actinides

"""

    strelements

Tuple containing all elements in order of increasing Z (up to Z=96) as a `String` followed by a space. Derived from [`SolarChem.periodictable`](@ref) for use inside `loadastromatdata`. A constant within the scope of the SolarChem module. 

"""
const strelements = string.(periodictable," ")


"""

    metadata = (:sample,:type,:group)

All metadata names in the NamedTuple generated from Astromat data. 
*Still in development* and eventually will become a constant in the SolarChem module.

"""
metadata = (:sample,:type,:group)


#=


"""
loadastromatdata(file::String)
"""



file = "../data/astromat/astromat-1_24_2024-eventId-1868.csv"

x = readdlm(string(@__DIR__,"/",file),',')


# headings = ("sample", "sample_url", "sample_url", "dataset_url", "dataset", "citation_url", "dataset_url", "citation", "citation_url", "collectionTypeName", "taxonType", "taxonName", "analysisType", "analyzedMaterialName", "analysisComment", "calcAvg", "numberOfReplicates")
headers = x[1,:]

sample = x[2:end,findfirst(x -> x=="sample", headers)]

taxonType = x[2:end,findfirst(x -> x=="taxonType", headers)]

taxonName = x[2:end,findfirst(x -> x=="taxonName", headers)]

analysisComment = x[2:end,findfirst(x -> x=="analysisComment", headers)]

citurl = x[2:end,findfirst(x -> x=="citation_url", headers)]



# analyzedMaterialName = x[2:end,findfirst(x -> x=="analyzedMaterialName", headers)] all "WHOLE ROCK"



# from_astromat = [:Fe, :As, :Cu, :Mn, :Ca, :Mg, :Zn, :S, :Al, :FeOT, :NiO, :FeS, :Fe2O3, :SO2, :Cr, :Ti, :Na, :P, :CARBONATE, :Fe2P, :Fe3P, :CoO, :K, :Cd, :Fe2O3T, :LOI, :Si, :Ag, :MoO3, :SO3, :V2O3, :BaO, :ST, :TiO, :ZrO2, :Ir, :Ga, :Ni, :Au, :Cr, :Co, :As, :Ge, :W, :Re, :Pt, :Cu, :Sb, :Sc, :K, :Zn, :V, :Se, :Os, :Hf, :Na, :Sr, :Ca, :C, :N, :Fe, :Mn, :Ru, :Rb, :Th, :U, :Zr, :Br, :Al, :Mg, :Ta, :Cs, :Pd, :Ti, :Te, :Y, :Cd, :Nb, :Mo, :Li, :Ag, :Tl, :Bi, :In, :Pb, :Sn, :P, :I, :S, :Si, :Be, :Rh, :B, :Hg, :Ba, :Sm, :Lu, :Yb, :La, :Eu, :Nd, :Cd, :Tb, :Ba, :Dy, :Ho, :Gd, :Er, :Pr, :Tm]

uncs= ["1S", "REL", "2S", "2S-ABS", "S-REL"]

oxides = [:CaO, :Na2O, :NiO, :CoO, :MgO, :FeO, :SiO2, :Al2O3, :MnO, :Cr2O3, :TiO2, :K2O, :P2O5]

elements = [:Li, :Be, :B, :C, :N, :Na, :Mg, :Al, :Si, :P, :S, :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Rb, :Sr, :Y, :Zr, :Nb, :Mo, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn, :Sb, :Te, :I, :Cs, :Ba, :La, :Pr, :Nd, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, :Hf, :Ta, :W, :Re, :Os, :Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Th, :U]

extras = [:FeOT, :FeS, :Fe2O3, :SO2, :CARBONATE, :Fe2P, :Fe3P, :Fe2O3T, :LOI,:MoO3, :SO3, :V2O3, :BaO, :ST, :TiO, :ZrO2, ]

volatiles = [:H, :CO2, :F, :H2OM, :H2O, :H2OP, :N2, :O, :HCl, :O2, :CO1  ]

exclude = ["lab", "method", "comment"] 

units = ["wt%", "%", "vol%", "mg/g", "μg/g", "ppm", "%m/m", "ppb" , "at%", "pg/g", "ng/g"]

volunites = ["umol/g", "ccstp", "ccstp/g"]


=#



