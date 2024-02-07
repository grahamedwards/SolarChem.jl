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

