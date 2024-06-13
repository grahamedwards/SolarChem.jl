"""
    @silence expr

Execute `expr` without printing any messages or warnings. If an error is encountered, it re-runs `expr` unsilenced.

# Examples

    julia> d = @silence sin(π) # place any variable assignment outside the macro.
    0.0

    julia> d
    0.0

    julia> @silence error("oh dear...")

    ┌ Warning: Error occurred within a @silence block. Rerunning unsilenced...
    └ @ Main ~/code/IsoMix.jl/test/silence.jl:21
    ERROR: oh dear...
    ...

"""
macro silence(block)
    quote
        ose,oso = stderr,stdout
        redirect_stdout(devnull)
        redirect_stderr(devnull)
        
        try 
            x = $(esc(block))
            redirect_stderr(ose)
            redirect_stdout(oso)
            x
        catch 
            redirect_stderr(ose)
            redirect_stdout(oso)
            @warn "Error occurred within a @silence block. Rerunning unsilenced..."
            $(esc(block))
        end
    end
end