# Auxilary functions
# Arc functions ────────────────────────────────────────────────────────────────────────────────
"""
    cₐ(a::Arc, x=a.x)

Returns arc cost for arc `a` for arc flow `x`
"""
function cₐ(a::Arc, x=a.x)
    tₒ= a.tₒ
    α = a.α
    β = a.β
    V = a.V
    ϕ = a.ϕ

    t = tₒ * (1 + α * (x/V) ^ β)
    t′= ϕ == 0 || β == 0 ? 0.0 : tₒ * α * β * (x ^ (β - 1))/(V ^ β)
    
    c = t + ϕ * x * t′
    return c
end

"""
    cₐ′(a::Arc, x=a.x)

Returns first derivative of arc cost wrt arc flow for arc `a` at arc flow `x`
"""
function cₐ′(a::Arc, x=a.x)
    tₒ= a.tₒ
    α = a.α
    β = a.β
    V = a.V
    ϕ = a.ϕ

    t′= β == 0 ? 0.0 : tₒ * α * β * (x ^ (β - 1))/(V ^ β)
    t″ = ϕ == 0 || β == 0 || β == 1 ?  0.0 : tₒ * α * β * (β - 1) * (x ^ (β - 2))/(V ^ β)
    
    c′ = t′ + ϕ * x * t″
    return c′
end



# Djikstra functions ────────────────────────────────────────────────────────────────────────────────
# NOTE: Node values must match node location (index) in the set of nodes N

"""
    djk(G::Graph, o::Origin)

Djikstra's label setting algorithm

Returns predecessor arc label (index in vector `A`) `Lᵖ` for least cost paths from origin `o` on graph `G`
"""
function djk(G::Graph, o::Origin)
    N = G.N
    A = G.A
    K = G.K
    r = o.n

    Lᵖ = Vector{Int64}(undef, length(N))     # Predecessor label
    Lᶜ = [Inf for _ in N]                    # Cost label
    Lᵖ[r.v] = -1
    Lᶜ[r.v] = 0.0
    
    X = copy(N)                              # Set of open nodes
    
    index = findfirst(x -> (x == r), X)::Int64
    t = X[index]
    deleteat!(X, index)
    while !isempty(X)
        i = t.v
        for j in t.H
            k = K[i,j]
            a = A[k]
            l = Lᶜ[i] + a.c
            if l ≤ Lᶜ[j] Lᵖ[j], Lᶜ[j] = k, l 
            end
        end
        index = argmin([Lᶜ[n.v] for n in X])
        t = X[index]
        deleteat!(X, index)
    end

    return Lᵖ
end

"""
    tree(G::Graph, o::Origin)

Returns tree for graph `G` rooted at origin `o`
"""
function tree(G::Graph, Lᵖ)
    N = G.N
    A = G.A

    T = Vector{Node}[[] for _ in N]

    for k in Lᵖ
        if k == -1 continue end
        a = A[k]
        t, h = a.t, a.h
        i = t.v
        push!(T[i], h)
    end

    return T
end

"""
    path(G::Graph, Lᵖ, r::Node, s::Node)

Returns path between origin node `r` and destination node `s` on graph `G` using predecessor label `Lᵖ` (for origin node `r`)
"""
function path(G::Graph, Lᵖ, r::Node, s::Node)
    A = G.A

    p = Arc[]

    h = s
    j = h.v
    while j ≠ r.v
        k = Lᵖ[j]
        a = A[k]
        push!(p, a)
        j = a.t.v
    end

    reverse!(p)
    return p
end
# ────────────────────────────────────────────────────────────────────────────────
