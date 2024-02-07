"""

    SolarChem.periodictable()

Return a Tuple containing all elements (as `Symbol`s) in order of increasing Z, up to Z=96.

"""
periodictable() = (
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

    SolarChem.strelements

A Tuple containing all elements in order of increasing Z (up to Z=96) as a `String` followed by a space. Derived from [`SolarChem.periodictable`](@ref) for use inside `loadastromatdata`. A constant within the scope of the SolarChem module. 

"""
const strelements = string.(periodictable()," ")


"""

    metadata = (:sample,:type,:group)

Return a Tuple of metadata names in the NamedTuple generated from Astromat data. 
*Still in development*.

"""
metadata() = (:sample,:type,:group)


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