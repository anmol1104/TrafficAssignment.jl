"""
    build(network, assignment)

Returns network as a graph with nodes, arcs, and relevant properties, for `:UE` or `:SO` assignment
"""
function build(network, assignment)
    #= ────────────────────────────────────────────────────────────────────────────────
    # NOTE: 
        Node values must match node location (index) in the set of nodes N. This 
        requires node values to start from 1 increasing linearly to the maximum value.
    ──────────────────────────────────────────────────────────────────────────────── =#

    ϕ = assignment == :UE ? false : true

    # network file
    ntwkfile = joinpath(@__DIR__, "Network\\$network\\network.csv")
    csv₁ = CSV.File(ntwkfile, types=[Int64, Int64, Float64, Float64, Float64, Float64, Float64])
    df₁ = DataFrame(csv₁)
    tail = df₁[!, 1]::Array{Int64,1}
    head = df₁[!, 2]::Array{Int64,1}
    linkcapacity = df₁[!, 3]::Array{Float64,1}
    linklength = df₁[!, 4]::Array{Float64,1}
    linkfft = df₁[!, 5]::Array{Float64,1}
    alpha = df₁[!, 6]::Array{Float64,1}
    beta = df₁[!, 7]::Array{Float64,1}

    # demand file
    dmndfile = joinpath(@__DIR__, "Network\\$network\\demand.csv")
    csv₂ = CSV.File(dmndfile)
    df₂ = DataFrame(csv₂)
    origin = df₂[!, 1]::Array{Int64,1}
    destination = df₂[!, 2]::Array{Int64,1}
    flow = df₂[!, 3]::Array{Float64,1}
    
    # Nodes
    N = Vector{Node}(undef, max(maximum(head), maximum(tail)))
    T = Vector{Int64}[[] for _ in 1:length(N)]
    H = Vector{Int64}[[] for _ in 1:length(N)]
    for row in 1:nrow(df₁)
        i = tail[row]
        j = head[row]
        push!(H[i], j)
        push!(T[j], i)    
    end
    for k in 1:length(N)
        n = Node(k, T[k], H[k])
        N[k] = n
    end

    # Arcs
    A = Vector{Arc}(undef, nrow(df₁))  
    xʳ= zeros(nrow(df₂))
    for row in 1:nrow(df₁)
        i, j = tail[row], head[row]
        t, h = N[i], N[j]
        V = linkcapacity[row]
        d = linklength[row]
        tₒ= linkfft[row]
        α = alpha[row]
        β = beta[row]
        a = Arc(t, h, V, d, tₒ, α, β, copy(xʳ), 0.0, 0.0, ϕ)
        A[row] = a
    end

    # Index mapping
    K = Dict((tail[row], head[row]) => row for row in 1:nrow(df₁))

    # Origins
    R = [N[k] for k in unique(origin)]
    Lᵖ= Vector{Int64}(undef, length(N))
    O = Vector{Origin}(undef, length(R))                            
    S = Vector{Node}[[] for _ in R]
    Q = Vector{Float64}[[] for _ in R]
    for row in 1:nrow(df₂)
        v = origin[row]
        r = N[v]
        v = destination[row]
        s = N[v]
        q = flow[row]
        k = findfirst(x -> (x == r), R)::Int64
        push!(S[k], s)
        push!(Q[k], q)
    end
    for k in 1:length(O)
        r = R[k]
        o = Origin(r, k, S[k], Q[k], copy(Lᵖ))
        O[k] = o
    end
    

    G = Graph(network, N, A, K, O)
    return G 
end