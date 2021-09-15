module TrafficAssignment

include("TAFW\\TAFW.jl")
include("TAPAS\\TAPAS.jl")

using .TAFW
using .TAPAS
using DataFrames
using Plots

"""
    assigntraffic(network; method=:FW, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:off)

Fukushima Frank-Wolfe method for traffic assignment.

# Returns
a named tuple with keys `:metadata`, `:report`, and `:output`
- `metadata::String`  : Text defining the traffic assignment run 
- `report::DataFrame` : A log of total network flow, total network cost, and run time for every iteration
- `output::DataFrame` : Flow and cost for every arc from the final iteration

# Arguments
- `network::String`         : Network
- `method::Symbol=:FW`      : One of `:FW`, `:fukushimaFW`, `:conjugateFW`
- `assignment::Symbol=:UE`  : Assignment type; one of `:UE`, `:SO`
- `tol::Float64=1e-5`       : Tolerance level for relative gap
- `maxiters::Int64=20`      : Maximum number of iterations
- `maxruntime::Int64=300`   : Maximum algorithm run time (seconds)
"""
function assigntraffic(network; method=:FW, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:off)
    if method == :FW  
        G = TAFW.build(network, assignment)
        return FW(G, tol, maxiters, maxruntime, log)
    elseif method == :fukushimaFW 
        G = TAFW.build(network, assignment)
        return fukushimaFW(G, tol, maxiters, maxruntime, log)
    elseif method == :conjugateFW 
        G = TAFW.build(network, assignment)
        return conjugateFW(G, tol, maxiters, maxruntime, log)
    elseif method == :TAPAS
        G = TAPAS.build(network, assignment)
        return itapas(G, tol, maxiters, maxruntime, log)
    else return (metadata="", report=DataFrame(), solution=DataFrame())
    end
end



"""
    compare(network; methods, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:on)

Compare assignment methods.
"""
function compare(network; methods, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:on)
    fig = plot()
    
    for method in methods
        _, report, = assigntraffic(network; method=method, assignment=assignment, tol=tol, maxiters=maxiters, maxruntime=maxruntime, log=log)
        y = report[!,:LOG₁₀RG]
        x = 0:length(y)-1
        plot!(x,y, label=String(method))
    end
    display(fig)
end



"""
    getrg(network, solution::DataFrame)

Returns relative gap for `network` given traffic assignment `solution`.
"""
function getrg(network, solution::DataFrame)
    G = build(network)
    N, A, O, K = G.N, G.A, G.O, G.K
    
    for row in 1:nrow(solution)
        i = solution[row, :FROM]::Int64
        j = solution[row, :TO]::Int64
        x = solution[row, :FLOW]
        c = solution[row, :COST]
        k = K[i,j]
        a = A[k]
        a.x = x
        a.c = c
    end

    num = 0.0
    for (i,o) in enumerate(O)
        r = o.n 
        Lᵖ = djk(G, o)
        for (j,s) in enumerate(o.S)
            qᵣₛ = o.Q[j]
            pᵣₛ = path(G, Lᵖ, r, s)
            num += qᵣₛ * cₚ(pᵣₛ)
        end 
    end

    den = 0.0
    for (k,a) in enumerate(A) den += a.x * a.c end

    rg = 1 - num/den
    return rg
end

export assigntraffic

end
#= ────────────────────────────────────────────────────────────────────────────────
# TODO:
1. Finish jobs in FW
2. Complete TAPAS
3. Complete tests
──────────────────────────────────────────────────────────────────────────────── =#