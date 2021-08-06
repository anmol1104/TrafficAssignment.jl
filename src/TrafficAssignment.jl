module TrafficAssignment

using CSV
using DataFrames
using Dates
using Printf
using Plots



struct Graph
    name::String                            # Network name
    N::Array{Int64,1}                       # Vector of nodes
    A::Array{Array{Int64,1},1}              # Vector of arcs as adjanceny list
    V::Array{Array{Float64,1},1}            # Vector of arc volume
    D::Array{Array{Float64,1},1}            # Vector of arc length
    T::Array{Array{Float64,1},1}            # Vector of arc free flow travel time
    α::Array{Array{Float64,1},1}            # Vector of arc bpr parameter
    β::Array{Array{Float64,1},1}            # Vector of arc bpr parameter
    R::Array{Int64,1}                       # Vector of origins
    S::Dict{Int64, Array{Int64,1}}          # Collection of destinations mapped to origin with non-zero flow (one to many mapping - r => {s; qᵣₛ ≠ 0})
    Q::Dict{Tuple{Int64, Int64}, Float64}   # Collection of flow mapped to origin-destination pair (one to one mapping - r,s => qᵣₛ; qᵣₛ ≠ 0)
end



include("build.jl")
include("segment.jl")
include("djk.jl")
include("FW.jl")
include("fukushimaFW.jl")
include("conjugateFW.jl")



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
    G = build(network)
    if method == :FW  return FW(G, assignment, tol, maxiters, maxruntime, log)
    elseif method == :fukushimaFW return fukushimaFW(G, assignment, tol, maxiters, maxruntime, log)
    elseif method == :conjugateFW return conjugateFW(G, assignment, tol, maxiters, maxruntime, log)
    else return (metadata="", report=DataFrame(), output=DataFrame())
    end
end



"""
    compare(network; methods, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:on)

Compare assignment methods.
"""
function compare(network; methods, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:on)
    fig = plot()
    
    G = build(network)
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
    N, A, R, S, Q = G.N, G.A, G.R, G.S, G.Q                                         # Graph
    x = [zeros(length(A[i])) for i in N]                                            # Arc flow
    c = [zeros(length(A[i])) for i in N]                                            # Arc cost
    
    for row in 1:nrow(solution)
        i = solution[row, :FROM]::Int64
        j = solution[row, :TO]::Int64
        k = findfirst(x -> (x == j), A[i])::Int64
        x[i][k] = solution[row, :FLOW]
        c[i][k] = solution[row, :COST]
    end

    num = 0.0
    for r in R 
        Lᵖ = djk(G, c, r)
        for s in S[r] 
            num += Q[r,s] * cₑ(G, c, path(Lᵖ, r, s))
        end 
    end

    den = 0.0
    for i in N
        for k in 1:length(A[i])
            den += x[i][k] * c[i][k]
        end 
    end

    rg = 1 - num/den
    return rg
end


export assigntraffic

end
#= ────────────────────────────────────────────────────────────────────────────────
# TODO
1. Complete tests
2. Check performance of conjugate fw againts pure fw and fukushima fw
3. Add tapas
──────────────────────────────────────────────────────────────────────────────── =#
