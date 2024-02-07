"""

    assigntype(str)

    Assign all corresponding genetic groups to a String `str` containing meteorite group/type information. See [`innergroups`](@ref) and [`outergroups`](@ref) for groups.

"""
function assigntype(v::String)
    
    i = findfirst(isdigit,v)

    if isnothing(i)
        i = 1
        v = ifelse(contains(v,"IMP"), "7", "0")
    end

    c = c2 = parse(Int,v[i])
()
    if contains(v, '-') 
        c2 = parse(Int,v[findfirst(x -> x=='-',v )+1])
    end

    if contains(v,'/')
        ic2 = findlast(x -> x=='/',v)+1
        if isdigit(v[ic2])
            c2 = parse(Int,v[ic2])
        end
    end   
    (c:c2...,)
end



"""

```julia
assigngroup(str)
```

Assign all corresponding genetic groups to a String `str` containing meteorite group/type information. See [`innergroups`](@ref) and [`outergroups`](@ref) for groups.

"""
function assigngroup(v::String)
    s = ()
    s = ifelse(isEH(v), (s..., :EH), s)
    s = ifelse(isEL(v), (s..., :EL), s)
    s = ifelse(isE(v), (s..., :E), s)
    s = ifelse(isH(v), (s..., :H), s)
    s = ifelse(isL(v), (s..., :L), s)
    s = ifelse(isLL(v), (s..., :LL), s)
    s = ifelse(isR(v), (s..., :R), s)
    s = ifelse(isCI(v), (s..., :CI), s)
    s = ifelse(isCV(v), (s..., :CV), s)
    s = ifelse(isCM(v), (s..., :CM), s)
    s = ifelse(isCR(v), (s..., :CR), s)
    s = ifelse(isCH(v), (s..., :CH), s)
    s = ifelse(isCB(v), (s..., :CB), s)
    s = ifelse(isCK(v), (s..., :CK), s)
    s = ifelse(isCL(v), (s..., :CL), s)
    s = ifelse(isCO(v), (s..., :CO), s)
    s = ifelse(isCung(v), (s..., :Cung), s)
    s
end

## Chondrite group is_ functions:

isEH(v::String) = contains(v, "EH")

isEL(v::String) = 0 < sum(contains.(v,("EL3", "EL3.", "EL4", "EL5", "EL6", "EL7", "EL/", "/EL", "EL(", "(EL)"))) 

function isE(v::String)
    c = ("E3", "E3.", "E4", "E5", "E6", "E7", "ENST")
    b = 0 < sum(contains.(v, c)) 
    b &= !contains(v,"ACHON") # exclude aubrites
end

function isH(v::String)
    c = ("H3", "H3.", "H4", "H5", "H6", "H7", "H/", "/H", "H-")
    b = 0 < sum(contains.(v, c)) 
    b *= !contains(v,"EH") 
    b *= !contains(v,"CH") #none in astromat when built but this could always change
end

function isL(v::String)
    c = ("L3", "L3.", "L4", "L5", "L6", "L7", "L/", "/L", "L(", "(L)")
    b = 0 < sum(contains.(v, c))  
    b *= !contains(v,"EL") 
    b *= !contains(v,"CL") #none in astromat when built but this could always change
    b *= ((contains(v,"/") | contains(v,")")) | !contains(v,"LL"))
end

isLL(v::String) = 0 < sum(contains.(v, ("LL3", "LL3.", "LL4", "LL5", "LL6", "LL7", "LL/", "/LL", "LL(", "(LL)")))

isR(v::String) = 0 < sum(contains.(v, ("R3", "R3.", "R4", "R5", "R6", "R7", "R/", "R(", "(R)")))

isCI(v::String) = contains(v, "CI1") 

isCV(v::String) = contains(v, "CV")

isCM(v::String) = contains(v, "CM") & !contains(v, "EUCRITE")

isCR(v::String) = contains(v, "CR") & !contains(v, "EUCRITE")

isCH(v::String) = 0 < sum(contains.(v, ("CH1", "CH2", "CH3")))

isCB(v::String) = contains(v, "CB")

isCK(v::String) = 0 < sum(contains.(v, ("CK3", "CK4", "CK5", "CK6", "CK7")))

isCL(v::String) = contains(v, "CL") & !contains(v, "UNCL")

isCO(v::String) =  0 < sum(contains.(v, ("CO1", "CO2", "CO3")))

isCung(v::String) = 0 < sum(contains.(v,("C-UNG", "C1", "C2", "C3", "C4")))

# unique(taxonName) =["LL3 CHONDRITE", "H4", "H3", "H6", "C-UNG", "H5", "L3", "CM", "DIOGENITE", "EUCRITE-MMICT", "EUCRITE", "H3.3 CHONDRITE", "CH3", "H3.0", "CK5", "EL3 CHONDRITE", "CM2", "ANGRITE", "EUCRITE-CM", "ACHONDRITE-UNG", "CR2", "", "UREILITE", "R4", "STONE-UNCL", "LUNAR", "EH4 CHONDRITE", "CV3", "CO3.5", "R3.8-5", "CO3.7", "L6", "E4", "EL5 CHONDRITE", "CI1", "L/LL4", "IRON-IVA", "LL5 CHONDRITE", "LL3.2/3.4 CHONDRITE", "EUCRITE (POLYMICT)", "H6 CHONDRITE", "L6 CHONDRITE", "L5 CHONDRITE", "CO3.6 CHONDRITE", "SHERGOTTITE", "L3.5 CHONDRITE", "L3.4 CHONDRITE", "L3.7-3.9 CHONDRITE", "MESOSIDERITE", "L3-6 CHONDRITE", "IRON-UNGROUPED", "EH3 CHONDRITE", "H3.7 CHONDRITE", "LL3.7 CHONDRITE", "CO3.0 CHONDRITE", "H3.9 CHONDRITE", "AUBRITE", "NIPR COLLECTION", "EUCRITE (UNBRECCIATED)", "CV3 CHONDRITE", "LUNAR-ANORTH. BRECCIA", "EL6 CHONDRITE", "ACAPULCOITE", "CM2 CHONDRITE", "CO3.4 CHONDRITE", "H5 CHONDRITE (IN ICE)", "H5 CHONDRITE", "UREILITE (AUG-BEARING)", "CK4 CHONDRITE", "L4 CHONDRITE", "LL3.2/3.5 CHONDRITE", "L3.4-3.7 CHONDRITE", "LL3.3 CHONDRITE", "LL6 CHONDRITE", "CM1/2 CHONDRITE", "CO3.5 CHONDRITE", "MARTIAN (OPX)", "SNC ORTHOPYROXENITE", "H4 CHONDRITE", "BRACHINITE", "L3.8 CHONDRITE", "LL3.4 CHONDRITE", "IRON-IIIAB", "L3.9 CHONDRITE", "EUCRITE (MG-RICH)", "H3.5-4 CHONDRITE", "CH3 CHONDRITE", "R3.6 CHONDRITE", "H4/5", "H3.5 CHONDRITE", "H3.4 CHONDRITE", "H5/6", "CHONDRITE-UNCL", "L3.7 CHONDRITE", "IRON, IAB-MG", "EUCRITE-PMICT", "L5", "CO3.6", "MARTIAN (SHERGOTTITE)", "L3.6 CHONDRITE", "MESOSIDERITE-B1", "L4", "CO3.0", "HOWARDITE", "L3.2 CHONDRITE", "LL3.5 CHONDRITE", "L3.7-4 CHONDRITE", "EUCRITE-UNBR", "CV3-AN", "LUNAR (ANORTH)", "H3.6 CHONDRITE", "LUNAR (GABBRO)", "EH4", "C2-UNG", "LL4 CHONDRITE", "EH4/5", "LL3.15", "LL3.9", "H/L3.9", "OC", "LL6 CHON. (BRECCIA)", "CR2 CHONDRITE", "CO3 CHONDRITE", "L/LL6", "CK3", "IRON-IVB", "CO3", "R3.8", "CK4", "MARTIAN (CHASSIGNITE)", "R3.5-6", "C3-UNG", "CK4/5", "EL4 CHONDRITE", "UREILITE-PMICT", "R3-4", "CK3-AN", "H-IMP MELT", "H3.8 CHONDRITE", "H~6", "LUNAR-BASALT", "H5-6 CHONDRITE", "CO3 CHONDRITE (ANOMALOUS)", "H5-7", "IRON, IIAB", "LUNAR (FELDSP. BRECCIA)", "UREILITE (POLYMICT)", "CK5 CHONDRITE", "IRON-IAB", "CM1 CHONDRITE", "C2 CHONDRITE UNGROUPED", "IRON-IIE (ANOMALOUS)", "L3.3-3.6 CHONDRITE", "ACAPULCOITE/LODRANITE", "IRON-IAB (ANOMALOUS)", "LUNAR-BASALTIC BRECCIA", "EUCRITE (BRECCIATED)", "UREILITE (POLYMICT ?)", "CK5/6 CHONDRITE", "L3.0 CHONDRITE", "CV3 CHONDRITE (REDUCED)", "H4 CHONDRITE (ANOMALOUS)", "EH4/5 CHONDRITE", "L/LL3.2 CHONDRITE", "TERRESTRIAL ROCK", "IRON-IIE", "CO3.3", "LODRANITE", "CR", "ACHON. UNGROUPED", "DIOGENITE (OLIVINE)", "EL4/5", "IRON-OCTAHEDRITE", "CHONDRITE UNGROUPED", "DIOGENITE (UNIQUE)", "CM CHONDRITE (ANOMALOUS)", "CR1 CHONDRITE", "LL3.8 CHONDRITE", "IRON, IAB-UNG", "L5-6", "C4-UNG", "WINONAITE", "CBB", "EL6/7", "L3.7-6", "R3-6", "EL7", "L/LL3.4", "CO3.8", "CO3.2", "K3", "L/LL5", "LL3.2", "R CHONDRITE", "H CHONDRITE (IMPACT MELT)", "L CHONDRITE (IMPACT MELT)", "R4 CHONDRITE", "CM1-2 CHONDRITE", "AUBRITE (ANOMALOUS)", "R6 CHONDRITE", "LUNAR (BASALT)", "LL5", "UREILITE (ANOMALOUS)", "CK6 CHONDRITE", "C3 CHONDRITE UNGROUPED", "L3.1 CHONDRITE", "L3.2-3.5 CHONDRITE", "E3 CHONDRITE (ANOMALOUS)", "CHON. (KAKANGARI-LIKE)", "EH5 CHONDRITE", "L5/6", "L(LL)3.05 CHONDRITE", "CV3 CHONDRITE (ANOMALOUS)", "R3 CHONDRITE", "L3.10 CHONDRITE", "NAKHLITE", "CB CHONDRITE", "LUNAR-FELDSPATHIC BRECCIA", "H5 CHONDRITE (ANOMALOUS)", "CR3 CHONDRITE", "CHONDRITE-UNG", "CM1", "CV2", "MARTIAN (NAKHLITE)", "IRON, IID-AN", "LL3.6 CHONDRITE", "CM-AN", "UNKNOWN", "MARTIAN (POLYMICT BRECCIA)", "R3/4", "R3", "CK5/6", "R3.9", "R5", "DIOGENITE-AN", "H", "ENST ACHON", "MARTIAN", "EUCRITE-AN", "DIOGENITE-PM", "LL5/6 CHONDRITE", "LL3.1", "LL3.05", "DIOGENITE-OLIVINE", "R3.7", "MARTIAN (AUGITE BASALT)", "H4-6 CHONDRITE", "C3", "CM1/2", "ACHONDRITE-PRIM", "E6", "CO3.4", "CH CHONDRITE", "CK3 CHONDRITE", "R3.8-6 CHONDRITE", "CBA", "ENSTATITE METEORITE UNGR", "CBB CHONDRITE", "CV3.4", "EUCRITE-BR", "H3.2-3.7", "MESOSIDERITE-AN", "L3.6-4 CHONDRITE", "R3.8-6", "R3-5", "KREEP BASALT", "H-METAL", "LL3.00", "L6/7", "H3-6", "C", "H/L3.6", "LL3.1-3.5 CHONDRITE", "H3-4", "LL7 CHONDRITE", "LODRANITE-AN", "H-AN", "LL3.0", "H7", "L4/5", "EH", "LL", "IRON, IIIAB", "IRON, IIE-AN", "H3/4", "E5-AN", "IRON, IAB COMPLEX", "EH6-AN", "H3.2-AN", "C1/2-UNG", "EH-IMP MELT", "CK6", "EH6 CHONDRITE", "PALLASITE, UNGROUPED", "R6"]
