# SolarChem plotting data

Import packages and load in serialized database file.

```julia (editor=true, logging=false, output=true)
using SolarChem, Statistics, CairoMakie
x = SolarChem.loadastromatdata();
```
```julia (editor=true, logging=false, output=true)
d = estimateuncertainty(trimextremes(x),5);

d = exclude(d, :comment, ("weather", "fusion", "CAI", "matrix", "chondrule")); # Exclude certain topics/comments: weather, fusion,  
d = excludeheated(d); 
```
## Outer Solar System

```julia (editor=true, logging=false, output=true)
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
countmeasurements(doss,allsolar())




```
```julia (editor=true, logging=false, output=true)
bse = bootstrapelements(100_000, doss, allsolar(Fe=true))
hh = histpanels(bse,labelsuffix=" (g/g)", darkmode=true)
hh
```
```julia (editor=true, logging=false, output=true)
ndoss = countmeasurements(doss, allsolar(Fe=true))
```
```julia (editor=true, logging=false, output=true)
bsr = bootstrapratios(100_000, doss, allsolar(Fe=true),:Fe, fractional=true)
hh = histpanels(bsr, labelsuffix="/Fe",darkmode=true)
#Label(hh[0,:],"elements",color=:white,fontsize=24)
hh
```
```julia (editor=true, logging=false, output=true)
outerFe=ratiosummary(bse, bsr, ratiocounts=countratios(doss, allsolar(),:Fe), minratios=3)
```
## Inner solar system

```julia (editor=true, logging=false, output=true)
diss = pullgroup(d, SolarChem.innergroups())

## Remove mistakes:

diss.Ca[findfirst(diss.Ca .> 0.9)] = .0092
diss.Cr[findfirst(diss.Cr .> 0.3)] = .0031

countmeasurements(diss,allsolar())
```
```julia (editor=true, logging=false, output=true)
bse = bootstrapelements(100_000, diss, allsolar())
hh = histpanels(bse,labelsuffix=" (g/g)", darkmode=true)
hh
```
```julia (editor=true, logging=true, output=true)
bsr = bootstrapratios(100_000, diss, allsolar(Fe=true),:Fe, fractional=true)
hh = histpanels(bsr, labelsuffix="/Fe",darkmode=true)
hh
```
```julia (editor=true, logging=false, output=true)
ndiss = countmeasurements(diss, majors())
innerFe = ratiosummary(bse, bsr)
lf = LoddersFegley1998()

for i = axes(innerFe[2],1)
    if isnan(innerFe[2][i,1])
    innerFe[2][i,:] .= lf.inner[innerFe[1][i]] .* (1, 0.1)
    end
end
hcat(innerFe[1], innerFe[2])
```
