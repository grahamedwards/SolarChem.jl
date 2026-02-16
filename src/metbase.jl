"""

    cleanmetbase(x)

Remove redundant MetBase data entries from a SolarChem dataset by replacing any redundant values in MetBase with `NaN`. Identifies MetBase entries through `MetBase` mentions under the `dataset` field and "redundant" values are identified using the `≈` operator. 

Because this method removes approximate redundancies, it does not remove inaccurate or redundant data that is obscured by typos MetBase or Astromat Synthesis. Simiarly, if a value published in an abstract (e.g. LPSC) is modified and later published in a peer-reviewed publication but the value changes, this method will retain both values. 

See also: [`cleanmetbase!`](@ref)

"""
function cleanmetbase(x::NamedTuple)
    y = deepcopy(x)
    cleanmetbase!(y)
    y
end

"""

    cleanmetbase!(x)

In place version of [`cleanmetbase`](@ref)

"""
function cleanmetbase!(x::NamedTuple)
    eachname = unique(x.name)
    allmetbase = findall(x -> occursin("MetBase", x), x.dataset)

    for i = eachindex(eachname)

        nameinds = findall(x -> x == eachname[i], x.name)

        metbase = intersect(nameinds,allmetbase) 
        nonmetbase = setdiff(nameinds,metbase)

        #fullcuts = onecuts = 0

        @inbounds for j in metbase

            els = ()
            @inbounds for k in periodictable()
                els = (els..., ifelse(isnan(x[k][j]), (), (k,))...)
            end

            @inbounds for l = nonmetbase
                if x.citation[j] == x.citation[l]
                    @inbounds for el in els
                        x[el][j] = NaN
                    end
                    #fullcuts +=1
                    # Note that citation intersects are fewer than actual `fullcuts` because of duplicate citations in both sets. 
                        #intersect(x.citation[metbase_intersect], x.citation[nonmetbase])
                    break
                else
                    @inbounds for el in els
                        same = x[el][j] ≈ x[el][l]
                        x[el][j] = ifelse(same, NaN, x[el][j])
                        #onecuts +=same
                    end
                end
            end
        end
    end
    x
end


"""

    removemetbase(x)

Remove all MetBase data entries from a SolarChem dataset by replacing any elemental data values with `NaN`. Identifies MetBase entries through `MetBase` mentions under the `dataset` field.

See also: [`removemetbase!`](@ref)

"""
function removemetbase(x::NamedTuple)
    y = deepcopy(x)
    removemetbase!(y)
    y 
end


"""

    removemetbase!(x)

In place version of [`removemetbase`](@ref)

"""
function removemetbase!(x::NamedTuple)
     allmetbase = findall(x -> occursin("MetBase", x), x.dataset)
     @inbounds for i = allmetbase
        @inbounds for k in periodictable()
            x[k][i] = NaN
        end
     end
end




