"""
build(; network)

Returns network as a graph of nodes, arcs, and relevant properties
"""
function build(network)
    # network file
    ntwkfile = joinpath(@__DIR__, "Network\\$network\\network.csv")
    csv₁ = CSV.File(ntwkfile, types=[Int64, Int64, Float64, Float64, Float64, Float64, Float64])
    df₁ = DataFrame(csv₁)
    head = df₁[!, 1]::Array{Int64,1}
    tail = df₁[!, 2]::Array{Int64,1}
    linkcapacity = df₁[!, 3]::Array{Float64,1}
    linklength = df₁[!, 4]::Array{Float64,1}
    linkfft = df₁[!, 5]::Array{Float64,1}
    alpha = df₁[!, 6]::Array{Float64,1}
    beta = df₁[!, 7]::Array{Float64,1}
    n = max(maximum(head), maximum(tail))
    
    N = [i for i in 1:n]                      # Nodes
    A = [Int64[] for _ in 1:n]                # Arcs as adjacency list
    V = [Float64[] for _ in 1:n]              # Link volume capcity
    D = [Float64[] for _ in 1:n]              # Link length
    T = [Float64[] for _ in 1:n]              # Link free flow travel time
    α = [Float64[] for _ in 1:n]              # BPR parameters
    β = [Float64[] for _ in 1:n]              # BPR parameters
    
    for i in 1:nrow(df₁)
        push!(A[head[i]], tail[i])
        push!(V[head[i]], linkcapacity[i])
        push!(D[head[i]], linklength[i])
        push!(T[head[i]], linkfft[i])
        push!(α[head[i]], alpha[i])
        push!(β[head[i]], beta[i])
    end
    
    # demand file
    dmndfile = joinpath(@__DIR__, "Network\\$network\\demand.csv")
    csv₂ = CSV.File(dmndfile)
    df₂ = DataFrame(csv₂)
    origin = df₂[!, 1]::Array{Int64,1}
    destination = df₂[!, 2]::Array{Int64,1}
    flows = df₂[!, 3]::Array{Float64,1}
    
    R = unique(origin)                                                        # Origins
    S = Dict(r => Int64[] for r in R)                                         # Destinations for every origin
    Q = Dict((origin[i], destination[i]) => flows[i] for i in 1:nrow(df₂))    # Demand between OD pairs
    
    for i in 1:nrow(df₂) push!(S[origin[i]], destination[i]) end
    

    G = Graph(network, N, A, V, D, T, α, β, R, S, Q)
    return G
end