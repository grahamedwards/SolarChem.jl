### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ ff549afa-84cc-11f0-3c3c-2b34c3e18e70
try 
	using SolarChem, PlutoUI
	import GLMakie as Mke;
catch
	import Pkg
	Pkg.add("GLMakie")
	Pkg.add("PlutoUI")
	Pkg.add(url="https://github.com/grahamedwards/SolarChem.jl", rev="main")
	using SolarChem, PlutoUI
end

# ╔═╡ 402f9139-c4a8-40bb-b1bb-1f5ee4403876
d0 = SolarChem.loadastromatdata()

# ╔═╡ 8ebfb59e-b265-4cbe-9c16-9cdb36d4879a
d1 = SolarChem.estimateuncertainty(SolarChem.trimextremes(d0),5);

# ╔═╡ fe92d954-e36a-44ed-8ddc-ddb261836222
d2 = SolarChem.exclude(d1, :comment, ("weather", "fusion")); # Exclude certain topics/comments: weather, fusion,

# ╔═╡ 8be31a53-b8b3-41ef-b32c-856d367bc43c
allmet = SolarChem.excludeheated(d2)

# ╔═╡ 27b6830c-45f7-417a-9053-567d3264030d
@bind grp Select([SolarChem.innergroups()..., SolarChem.outergroups()...])

# ╔═╡ 866cf5c5-fb81-47fd-a982-16c840be5e8d
groupdata = pullgroup(allmet,grp)

# ╔═╡ a653c3b5-c9e5-42d9-8014-2ccd38d8ff9f
@bind el Select([SolarChem.periodictable()...])

# ╔═╡ 83257936-e505-4f39-bed7-bfb681d5581b
d = let
	d = groupdata[el]
	d[.!isnan.(d)]
end

# ╔═╡ 1c2ea3b1-0029-4719-aa7b-9a4a657c0e9b
@bind xmin PlutoUI.Slider(isempty(d) ? [0,1] : LinRange(minimum(d), maximum(d), 100))

# ╔═╡ a13bd2dd-d5c6-4dd7-8f83-a37301cd2125
@bind xmax PlutoUI.Slider(isempty(d) ? [0,1] : LinRange(xmin, maximum(d), 100), default=maximum(d))

# ╔═╡ 13c6169f-1643-48a6-a296-be6116ed4272
let 
	f = Mke.Figure()
	ax = Mke.Axis(f[1,1], limits=(xmin, xmax, nothing, nothing))

	if !isempty(d)
		Mke.hist!(ax,d)
	else
		Mke.text!(ax,"No data")
	end
	f
end

# ╔═╡ ca522569-8c91-4608-9a94-9879451b7c08
let 
	report = "| Name  | Citation |\n|:-- |:-- |\n"
	citations = ()
	for i = eachindex(groupdata.name)
		if xmin <= groupdata[el][i] <= xmax
			citations = (citations..., groupdata.citation[i])
			report *= string(groupdata.name[i], 
							 " | [", 
							 groupdata.citation[i], 
							 "](", 
							 groupdata.citation[i], ")\n")
		end
	end
	report *= "\n\nUnique citations\n\n"
	uc = unique(citations)
	 for i in eachindex(uc)
		 uci = uc[i]
		 report *= "[$uci]($uci)\n"
	 end
	Markdown.parse(report)
end

# ╔═╡ Cell order:
# ╠═ff549afa-84cc-11f0-3c3c-2b34c3e18e70
# ╠═402f9139-c4a8-40bb-b1bb-1f5ee4403876
# ╠═8ebfb59e-b265-4cbe-9c16-9cdb36d4879a
# ╠═fe92d954-e36a-44ed-8ddc-ddb261836222
# ╠═8be31a53-b8b3-41ef-b32c-856d367bc43c
# ╠═866cf5c5-fb81-47fd-a982-16c840be5e8d
# ╠═83257936-e505-4f39-bed7-bfb681d5581b
# ╟─1c2ea3b1-0029-4719-aa7b-9a4a657c0e9b
# ╟─a13bd2dd-d5c6-4dd7-8f83-a37301cd2125
# ╟─13c6169f-1643-48a6-a296-be6116ed4272
# ╟─27b6830c-45f7-417a-9053-567d3264030d
# ╟─a653c3b5-c9e5-42d9-8014-2ccd38d8ff9f
# ╟─ca522569-8c91-4608-9a94-9879451b7c08
