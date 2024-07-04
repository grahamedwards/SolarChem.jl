## Make NamedTuples look nice. Should be fairly harmless type piracy :)
function Base.display(x::NamedTuple)
    println("NamedTuple with $(length(x)) elements:")
    kx = keys(x)
    @inbounds for i in eachindex(ifelse(length(kx) > 10, 1:10, kx))
        xi = x[i]
        print("  $(kx[i])")
        if xi isa NamedTuple
            println("::NamedTuple{$(eltype(xi))} of length $(length(xi))")
        elseif xi isa Composition
            println(" = Composition($(xi.m), $(xi.s))")
        elseif xi isa AbstractRange 
            println("::$(typeof(xi)) = $xi")
        elseif length(xi) > 5
            op, cl = xi isa Tuple ? ("(",")") : ("[","]")
            println("::$(typeof(xi)) = $op ", first(xi), "...", last(xi),cl)
        else
            println("::$(typeof(xi)) = $xi")
        end
    end
    length(x) > 10 && println("  ⋮")
end