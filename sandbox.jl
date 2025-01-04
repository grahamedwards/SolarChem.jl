# Sandbox

using SolarChem
import Serialization

x = loadastromatdata();
lf = LoddersFegley1998()
# Serialization.serialize("workshop/astromat.cereal", loadastromatdata())
# x = Serialization.deserialize("workshop/astromat.cereal")

d = estimateuncertainty(trimextremes(x),5);

d = exclude(d, :comment, ("weather", "fusion", "CAI", "matrix", "chondrule")); # Exclude certain topics/comments: weather, fusion,  
d = excludeheated(d); 



## Inner solar system:
diss = pullgroup(d, SolarChem.innergroups())
## Remove mistakes:
    diss.Ca[findfirst(diss.Ca .> 0.9)] = .0092
    diss.Cr[findfirst(diss.Cr .> 0.3)] = .0031
bse = bootstrapelements(100_000, diss, allsolar())
bsr = bootstrapratios(100_000, diss, allsolar(Fe=true),:Fe, fractional=true)
innerFe = ratiosummary(bse, bsr)


## Outer solar system
doss = pullgroup(d, SolarChem.outergroups())
## Remove outliers, errors:
    doss.C[findfirst(x -> x =="https://page.astromat.org/dataset/fd3cynsBjW5SsTLVHRpM", doss.dataset)]=NaN # Anomalous, high-C outlier.
    allendecai = findfirst(doss.name .== "Allende" .&& doss.dataset .== "https://page.astromat.org/dataset/T-yOt3QBpybrwNsBY-ax")
    anomhighTi = findfirst(doss.dataset.=="https://page.astromat.org/dataset/ktygynsBjW5SsTLVdPyj") # Anomolous chondrite ALH85085 with very high bulk Ti. 

    pcp = findall(doss.dataset .== "https://page.astromat.org/dataset/HwqTWoYB_AwVtctyi_xc") # "Poorly characterized phases" or serpentine-tochilinite intergrowths.

    for i = keys(doss)
        if i ∉ SolarChem.metadata()
            doss[i][373]= doss[i][1503]=  doss[i][findfirst(doss.dataset .== "https://page.astromat.org/dataset/Uw23aYYB_AwVtctyZ9Dk")] = NaN # AstroMat errors.
            doss[i][anomhighTi] = doss[i][allendecai] = NaN 
            doss[i][pcp] .= NaN
        end
    end

bse = bootstrapelements(100_000, doss, allsolar(Fe=true))
bsr = bootstrapratios(100_000, doss, allsolar(Fe=true),:Fe, fractional=true)

outerFe = ratiosummary(bse, bsr)


## Run the metropolis

iss, oss = NamedTuple(innerFe), NamedTuple(outerFe)
iss, oss = removefrom(:C, iss), removefrom(:C, iss);





p = Fractions(.2,.2)
j = Fractions(.001,0.001)
out = solarmixmetropolis(10_000, p, j, iss, oss, burnin=10_000)