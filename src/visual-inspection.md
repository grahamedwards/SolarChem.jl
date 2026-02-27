# SolarChem visual inspection

Welcome! This notebook provides an interactive interface for examining data pulled from NASA's Astromaterials Database, [Astromat](https://astromat.org/).

## Getting started.

Click the play button above to run all the cells in the notebook. This will get everything loaded and may take a couple minutes, especially the first time you use this notebook. 

To run a cell, press `Shift+Enter`. You may find it helps to re-run a cell to improve performance or update information.

### A word of caution

Both [SolarChem](https://github.com/grahamedwards/SolarChem.jl) (which provides the underlying maths at work here) and [BonitoBook](https://github.com/SimonDanisch/BonitoBook.jl/) (which provides this notebook environment) are new and experimental softwares still under active development. That means this notebook is a little rough around the edges. Please be patient as you encounter challenges, and I encourage you to read comments below (indicated with a `#`) to understand what is going on. Please reach out to me if you have problems or questions!

```julia (editor=true, logging=false, output=true)
try 
    import WGLMakie as Mke
catch 
    import Pkg
    Pkg.add("WGLMakie")
    import WGLMakie as Mke
end
using SolarChem
allmet = let
	d0 = fastromat()
	d1 = SolarChem.estimateuncertainty(SolarChem.trimextremes(d0),5)
	d2 = SolarChem.exclude(d1, :comment, ("weather", "fusion")); # Exclude certain topics/comments: weather, fusion,
    d3 = SolarChem.excludeheated(d2) # exclude comments that indicate heating experiments
	cleanmetbase(d3) # remove redundant metbase entries
end

elements = collect(periodictable()) # create a vector of elements from the `periodictable` function
groups = [SolarChem.innergroups()..., SolarChem.outergroups()...] # combine lists of inner and outer meteorite groups into a single vector.

"""
    printcitations(printcitations(x, ii, y, el)

    Print citations from a SolarChem data NamedTuple, `x`, a list of indices `ii`, a specific element of interest `el`, and the range of values for that element `y`.
"""
function printcitations(x, ii, y, el)
    sldr = y[1] .< x[el] .< y[2]
    i = ii .& sldr
    report = "| Name  | Citation |\n|:-- |:-- |\n"
    citations = ()
    nm, ct = x.name[i], x.citation[i]
    for i = eachindex(nm)
    citations = (citations..., ct[i])
    report *= string(nm[i], " | [`", ct[i], "`](", ct[i], ")\n")   
    end
    report *= "\n\nUnique citations\n\n"
    uc = unique(citations)
    for i in eachindex(uc)
        uci = uc[i]
        report *= "[`$uci`]($uci)\n"
    end
    Markdown.parse(report)
end
```
## Visual inspection

Run the cell below to generate the plot and then scroll down to check it out...

```julia (editor=true, logging=false, output=true)
f = Mke.Figure()
Mke.DataInspector(f)
ax = Mke.Axis(f[1, 1])
Mke.deactivate_interaction!(ax, :scrollzoom)
Mke.deactivate_interaction!(ax, :rectanglezoom)
domainslider = Mke.IntervalSlider(f[2, 1])

grpmenu = Mke.Menu(f, options = groups, direction = :up, default="LL")
elmenu = Mke.Menu(f, options=elements, direction=:down, default="Na")

f[1,2] = Mke.vgrid!(
    Mke.Label(f, "Group", width = nothing), grpmenu,
    Mke.Label(f, "Element", width = nothing), elmenu,
    tellheight = false, width = 100)

grp , el = Observable{Symbol}(), Observable{Symbol}()

xgrp, igl, xgl, igld, xgld = Observable{Any}(), Observable{Any}(), Observable{Any}(), Observable{Any}(), Observable{Any}()

Mke.on(grpmenu.selection) do g
    grp[] = Symbol(g)
end
Mke.notify(grpmenu.selection)

Mke.on(elmenu.selection) do elsec
    el[] = Symbol(elsec)
end
Mke.notify(elmenu.selection)

xgrp = Mke.lift(grp) do g 
    pullgroup(allmet,g)
end

igl = Mke.lift(xgrp, el) do x,l
    d = getproperty(x,l)
    .!isnan.(d) # indices corresponding to element and group. 
end

xgl = Mke.lift(xgrp, el, igl) do x,l,i
    xx = getproperty(x,l)[i] 
    Mke.autolimits!(ax)
    xmin, xmax = extrema(xx)
    domainslider.range[] = LinRange(xmin, xmax, 200)
    domainslider.interval[] = (xmin, xmax)
    xx
end

xgld = Mke.lift((x,y) -> x[y[1] .< x .< y[2]], xgl, domainslider.interval)

Mke.on(domainslider.interval) do x 
    Mke.ylims!(nothing, nothing)
end
Mke.notify(domainslider.interval)

Mke.hist!(xgld, color=:gray50)
f
```
## Instructions

The histogram above is generated using the built-in interactivity of the [Makie](https://docs.makie.org/stable/) plotting package.


### Group and element dropdowns
Use the dropdown menus on the right to select the meteorite group and element you wish to examine. 


### Domain slider
Use the slider below the x-axis to adjust the domain of your histogram.  Click the circles at either end and drag them to focus on a specific region.  Note that bins outside the domain disappear!

### Datatips 
Hovering over a histogram bar reveals a datatip that reports that bin's midpoint value (x) and count (y). While this does not provide the exact values within a bin, if you zoom enough for single count bins, it should provide a pretty close estimate.

### Troubleshooting

If the plot starts behaving strangely, try the following sequence:

1. Select a new group and a new element.
2. If that doesn't work, re-run the cell above (`Shift+Enter`) to refresh the plot. Unfortunately, this will reset the group and element to the default.
3. If the dropdown selectors are unable to get you to the element/group pairing of interest, you can manually reset the defaults in the code for the plot above (lines 8-9). If you have to resort to this, please let me know!

### Tabulated data

You can print the data shown on the plot by running the cell below (`Shift+Enter`). Remember to re-run the cell below after any interaction with the plot above to get the latest list of citations corresponding to visible data. 
```julia (editor=true, logging=false, output=true)
printcitations(xgrp[], igl[], domainslider.interval[], el[])
```
