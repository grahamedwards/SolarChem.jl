"""

    SolarChem.periodictable()

Return a Tuple containing all elements (as `Symbol`s) in order of increasing Z, up to Z=96.

"""
periodictable() = (
# row 1
    :H, :He, 
# row 2
    :Li, :Be, :B, :C, :N, :O, :F, :Ne, 
# row 3    
    :Na, :Mg, :Al, :Si, :P, :S, :Cl, :Ar, 
# row 4
    :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Kr, 
# row 5
    :Rb, :Sr, :Y, :Zr, :Nb, :Mo, :Tc, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn, :Sb, :Te, :I, :Xe,
# row 6
    :Cs, :Ba, # groups 1-2
    :La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu, :Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, # REEs/group 3
    :Hf, :Ta, :W, :Re, :Os, :Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Po, :At, :Rn, # groups 4-18
#row 7, some actinides
    :Fr, :Ra, :Ac, :Th, :Pa, :U, :Np, :Pu, :Am, :Cm) 

"""

    SolarChem.strelements

A Tuple containing all elements in order of increasing Z (up to Z=96) as a `String` followed by a space. Derived from [`SolarChem.periodictable`](@ref) for use inside `loadastromatdata`. A `const` within the scope of the SolarChem module. 

"""
const strelements = string.(periodictable()," ")



"""

    oxideconversion

A NamedTuple with names of elements corresponding to a multiplicative factor to convert oxide composition to elemental composition. `Fe` corresponds to ferrous iron (FeO) and `Fe3` to ferric iron (Fe₂O₃). `oxideconversion` is a `const` within the scope of the SolarChem module. 

---

    CaO * oxideconversion.Ca = Ca

    Ca / oxideconversion.Ca = CaO

"""
const oxideconversion = (; 
    Ca= 40.078/(40.078 + 15.999), 
    Na = 2* 22.989769 / (2* 22.989769 + 15.999), 
    Ni = 58.6934 / (58.6934 + 15.999), 
    Mg = 24.305 / (24.305 + 15.999), 
    Co = 58.933195 / (58.933195 + 15.999), 
    Fe = 55.845 / (55.845 + 15.999), 
    Si = 28.0855 / (28.0855 + 2* 15.999), 
    Al = 2* 26.981539 / (2* 26.981539 + 3* 15.999), 
    Mn = 54.938044 / (54.938044 + 15.999), 
    Cr = 2* 51.9961 / (2* 51.9961 + 3* 15.999), 
    Ti = 47.867 / (47.867 + 2* 15.999), 
    K = 2* 39.0983 / (2* 39.0983 + 15.999), 
    P = 2* 30.973762 / (2* 30.973762 + 5* 15.999), 
    Fe3 = 2* 55.845 / (2* 55.845 + 3* 15.999)
    )


"""

    metadata = (:name,:type,:group)

Return a Tuple of metadata names in the NamedTuple generated from Astromat data. 
*Still in development*.

"""
metadata() = (:name, :type, :group, :comment, :citation, :dataset)


"""

    innergroups()

Return a Tuple containing all chondrite genetic groups associated with the inner solar system and used in SolarChem.

"""
innergroups() = (:EH, :EL, :E, :H, :L, :LL, :R)

"""

    outergroups()

Return a Tuple containing all chondrite genetic groups associated with the outer solar system and used in SolarChem.

"""
outergroups() = (:CI, :CV, :CM, :CR, :CH, :CB, :CK, :CL, :CO, :Cung)



"""

    majors()

Return a Tuple containing major elements.

"""
majors() = (:Na, :Mg, :K, :Ca, :Al, :Fe, :Ni, :Si)



"""

    allsolar(; Fe=true)

Return a Tuple of all elements with solar twin data reported in Bedell+ 2018 (*ApJ*, [doi:10.3847/1538-4357/aad908](https://doi.org/10.3847/1538-4357/aad908)). Excludes `:Fe` by default, but `Fe=true` will include it. 

see also: [`lightsolar`](@ref), [`mediumsolar`](@ref), [`heavysolar`](@ref)

"""
allsolar(; Fe=true) = (lightsolar()...,mediumsolar(Fe=Fe)...,heavysolar()...)

"""

    lightsolar()

Return a Tuple of Period 1-3 elements with solar twin data reported in Bedell+ 2018 (*ApJ*, [doi:10.3847/1538-4357/aad908](https://doi.org/10.3847/1538-4357/aad908)).

see also: [`allsolar`](@ref)

"""
lightsolar() = (:C, :Na, :Mg, :Al, :Si, :S)

"""

    mediumsolar(; Fe=true)

Return a Tuple of Period 4 elements with solar twin data reported in Bedell+ 2018 (*ApJ*, [doi:10.3847/1538-4357/aad908](https://doi.org/10.3847/1538-4357/aad908)).

Excludes `:Fe` by default, but `Fe=true` will include it. 

see also: [`allsolar`](@ref)

"""
mediumsolar(;Fe=true) = ifelse(Fe, (:Ca, :Sc, :Ti, :V, :Cr, :Mn, :Fe, :Co, :Ni, :Cu, :Zn), (:Ca, :Sc, :Ti, :V, :Cr, :Mn, :Co, :Ni, :Cu, :Zn) )

"""

    heavysolar()

Return a Tuple of Period 5-6 elements with solar twin data reported in Bedell+ 2018 (*ApJ*, [doi:10.3847/1538-4357/aad908](https://doi.org/10.3847/1538-4357/aad908)).

see also: [`allsolar`](@ref)
"""
heavysolar() = (:Sr, :Y, :Zr, :Ba, :La, :Ce, :Pr, :Nd, :Sm, :Eu, :Gd, :Dy)

