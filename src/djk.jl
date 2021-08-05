"""
    djk(G, c, r)

Djikstra's label setting algorithm

Returns predecessor label `Lᵖ` for least cost paths from origin node `r` on graph `G` with cost structure `c`
"""
function djk(G::Graph, c, r)
    N, A = G.N, G.A
    Lᵖ = [if i == r r else -1 end for i in N]           # Predecessor label
    Lᶜ = [if i == r 0.0 else Inf end for i in N]        # Cost label
    X = copy(N)                                         # Set of open nodes
    
    i = r
    deleteat!(X, i)
    while !isempty(X)
        for (k,j) in enumerate(A[i])
            l = Lᶜ[i] + c[i][k]
            if l < Lᶜ[j] && j ∈ X 
                Lᵖ[j], Lᶜ[j] = i, l 
            end
        end
        index = argmin([Lᶜ[i] for i in X])
        i = X[index]
        deleteat!(X, index)
    end
    return Lᵖ
end

"""
    tree(Lᵖ)

Returns tree for graph `G` rooted at node `r` given predecessor label `Lᵖ` (for node `r`)
"""
function tree(G::Graph, Lᵖ)
    N = G.N
    T = Array{Int64,1}[[] for _ in N]
    for j in N
        i = Lᵖ[j]
        if i ≠ j && i ≠ -1 push!(T[i], j) end
    end
    return T
end

"""
    path(Lᵖ, r, s)

Returns path between origin node `r` and destination node `s` using predecessor label `Lᵖ` (for node `r`)
"""
function path(Lᵖ, r, s)
    p = Int64[]
    i = s
    push!(p, i)
    while i ≠ r
        i = Int(Lᵖ[i])
        push!(p, i)
    end
    reverse!(p)
    return p
end