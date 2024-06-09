using SolarChem, Test
using StableRNGs

# Note that silent errors happen within a @silence block.
macro silence(block)
    quote
        ose,oso = stderr,stdout
        redirect_stderr(devnull); redirect_stdout(devnull)
        x= $(esc(block))
        redirect_stderr(ose); redirect_stdout(oso)
        x
    end
end

mean(x) = SolarChem.StatsBase.mean(x)

@testset "Load Data" begin include("load-data.jl") end
@testset "Data Mgmt" begin include("data-mgmt.jl") end
@testset "Resamplers" begin include("resamplers.jl") end

