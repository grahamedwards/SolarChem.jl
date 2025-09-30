# Metropolis sandbox

using SolarChem
import Serialization

x = SolarChem.loadastromatdata();
d = estimateuncertainty(trimextremes(x),5);

d = exclude(d, :comment, ("weather", "fusion", "CAI", "matrix", "chondrule")); # Exclude certain topics/comments: weather, fusion,  
d = excludeheated(d); 

## Inner solar system:
diss = pullgroup(d, SolarChem.innergroups())
bse = bootstrapelements(100_000, diss, allsolar())
bsr = bootstrapratios(100_000, diss, allsolar(Fe=true),:Fe, fractional=true)
innerFe = ratiosummary(bse, bsr)

## Outer solar system
doss = pullgroup(d, SolarChem.outergroups())
bse = bootstrapelements(100_000, doss, allsolar(Fe=true))
bsr = bootstrapratios(100_000, doss, allsolar(Fe=true),:Fe, fractional=true)
outerFe = ratiosummary(bse, bsr)

## Run the metropolis
iss, oss = NamedTuple(innerFe), NamedTuple(outerFe)
iss, oss = removefrom(:C, iss), removefrom(:C, iss);

p = Fractions(.2,.2)
j = Fractions(.001,0.001)
out = solarmixmetropolis(10_000, p, j, iss, oss, burnin=10_000)