
"""
loadastromatdata(file::String)
"""

using DelimitedFiles, VectorizedStatistics


file = "../data/astromat/astromat-1_24_2024-eventId-1868.csv"

x = readdlm(string(@__DIR__,"/",file),',')

# from_astromat = [:Fe, :As, :Cu, :Mn, :Ca, :Mg, :Zn, :S, :Al, :FeOT, :NiO, :FeS, :Fe2O3, :SO2, :Cr, :Ti, :Na, :P, :CARBONATE, :Fe2P, :Fe3P, :CoO, :K, :Cd, :Fe2O3T, :LOI, :Si, :Ag, :MoO3, :SO3, :V2O3, :BaO, :ST, :TiO, :ZrO2, :Ir, :Ga, :Ni, :Au, :Cr, :Co, :As, :Ge, :W, :Re, :Pt, :Cu, :Sb, :Sc, :K, :Zn, :V, :Se, :Os, :Hf, :Na, :Sr, :Ca, :C, :N, :Fe, :Mn, :Ru, :Rb, :Th, :U, :Zr, :Br, :Al, :Mg, :Ta, :Cs, :Pd, :Ti, :Te, :Y, :Cd, :Nb, :Mo, :Li, :Ag, :Tl, :Bi, :In, :Pb, :Sn, :P, :I, :S, :Si, :Be, :Rh, :B, :Hg, :Ba, :Sm, :Lu, :Yb, :La, :Eu, :Nd, :Cd, :Tb, :Ba, :Dy, :Ho, :Gd, :Er, :Pr, :Tm]

uncs= ["1S", "REL", "2S", "2S-ABS", "S-REL"]

oxides = [:CaO, :Na2O, :NiO, :CoO, :MgO, :FeO, :SiO2, :Al2O3, :MnO, :Cr2O3, :TiO2, :K2O, :P2O5]

elements = [:Li, :Be, :B, :C, :N, :Na, :Mg, :Al, :Si, :P, :S, :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Rb, :Sr, :Y, :Zr, :Nb, :Mo, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn, :Sb, :Te, :I, :Cs, :Ba, :La, :Pr, :Nd, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, :Hf, :Ta, :W, :Re, :Os, :Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Th, :U]

extras = [:FeOT, :FeS, :Fe2O3, :SO2, :CARBONATE, :Fe2P, :Fe3P, :Fe2O3T, :LOI,:MoO3, :SO3, :V2O3, :BaO, :ST, :TiO, :ZrO2, ]

volatiles = [:H, :CO2, :F, :H2OM, :H2O, :H2OP, :N2, :O, :HCl, :O2, :CO1  ]

exclude = ["lab", "method", "comment"] 

units = ["wt%", "%", "vol%", "mg/g", "μg/g", "ppm", "%m/m", "ppb" , "at%", "pg/g", "ng/g"]

volunites = ["umol/g", "ccstp", "ccstp/g"]

periodictable = [
    :H, :He, 
    :Li, :Be, :B, :C, :N, :O, :F, :Ne, 
    :Na, :Mg, :Al, :Si, :P, :S, :Cl, :Ar, 
    :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Kr, 
    :Rb, :Sr, :Y, :Zr, :Nb, :Mo, :Tc, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn, :Sb, :Te, :I, :Xe, 
    :Cs, :Ba, :La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, :Hf, :Ta, :W, :Re, :Os, :Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Po, :At, :Rn, 
    :Fr, :Ra, :Ac, :Th, :Pa, :U, :Np, :Pu, :Am, :Cm]

    ["H ", "He ", "Li ", "Be ", "B ", "C ", "N ", "O ", "F ", "Ne ", "Na ", "Mg ", "Al ", "Si ", "P ", "S ", "Cl ", "Ar ", "K ", "Ca ", "Sc ", "Ti ", "V ", "Cr ", "Mn ", "Fe ", "Co ", "Ni ", "Cu ", "Zn ", "Ga ", "Ge ", "As ", "Se ", "Br ", "Kr ", "Rb ", "Sr ", "Y ", "Zr ", "Nb ", "Mo ", "Tc ", "Ru ", "Rh ", "Pd ", "Ag ", "Cd ", "In ", "Sn ", "Sb ", "Te ", "I ", "Xe ", "Cs ", "Ba ", "La ", "Ce ", "Pr ", "Nd ", "Pm ", "Sm ", "Eu ", "Gd ", "Tb ", "Dy ", "Ho ", "Er ", "Tm ", "Yb ", "Lu ", "Hf ", "Ta ", "W ", "Re ", "Os ", "Ir ", "Pt ", "Au ", "Hg ", "Tl ", "Pb ", "Bi ", "Po ", "At ", "Rn ", "Fr ", "Ra ", "Ac ", "Th ", "Pa ", "U ", "Np ", "Pu ", "Am ", "Cm "]




