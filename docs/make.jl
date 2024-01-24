using SolarChem
using Documenter

DocMeta.setdocmeta!(SolarChem, :DocTestSetup, :(using SolarChem); recursive=true)

makedocs(;
    modules=[SolarChem],
    authors="Graham Harper Edwards",
    sitename="SolarChem.jl",
    format=Documenter.HTML(;
        canonical="https://grahamedwards.github.io/SolarChem.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/grahamedwards/SolarChem.jl",
    devbranch="main",
)
