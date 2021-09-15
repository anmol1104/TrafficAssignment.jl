# Auxilary functions
# Arc functions ────────────────────────────────────────────────────────────────────────────────
"""
    cₐ(a::Arc, x=a.x)

Returns arc cost for arc `a` for arc flow `x`
"""
function cₐ(a::Arc, x=sum(a.xʳ))
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
function cₐ′(a::Arc, x=sum(a.xʳ))
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


# Segment functions ────────────────────────────────────────────────────────────────────────────────
"""
    cₑ(e::Vector{Arc})

Returns segment cost for segment `e`
"""
function cₑ(e::Vector{Arc})
    c = 0.0
    for a in e c += a.c end
    return c
end

"""
    cₑ′(e::Vector{Arc})

Returns first derivative of segment cost for segment `e`
"""
function cₑ′(e::Vector{Arc})
    c′ = 0.0
    for a in e c′ += a.c′ end
    return c′
end

"""
    fₑ(e::Vector{Arc}, o::Origin)

Return minimum flow on segment `e` from origin `o`
"""
function fₑ(e::Vector{Arc}, o::Origin)
    f = Inf
    k = o.k
    for a in e if a.xʳ[k] < f f = a.xʳ[k] end end
    return f
end


# Djikstra functions ────────────────────────────────────────────────────────────────────────────────
# NOTE: Node values must match node location (index) in the set of nodes N

"""
    djk(G::Graph, o::Origin)

Djikstra's label setting algorithm.
Returns predecessor arc label (index in vector `A`) `Lᵖ` for least cost paths from origin `o` on graph `G`.
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
            if l < Lᶜ[j] && N[j] ∈ X Lᵖ[j], Lᶜ[j] = k, l 
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
path(G::Graph, o::Origin, s::Node) = path(G, o.Lᵖ, o.n, s)  # Shorter version

# ────────────────────────────────────────────────────────────────────────────────

"""
    ispotential(a::Arc, o::Origin, G::Graph)

Identfies if arc `a` on graph `G` is a potential arc wrt flow from origin `o`
"""
function ispotential(a::Arc, o::Origin, G::Graph)
    k = o.k
    xʳₐ = a.xʳ[k]
    cₜₕ = a.c

    t = a.t
    pᵣₜ = path(G, o, t)
    uʳₜ = 0.0
    for a in pᵣₜ uʳₜ += a.c end
    
    h = a.h
    pᵣₕ = path(G, o, h)
    uʳₕ = 0.0
    for a in pᵣₕ uʳₕ += a.c end
    
    πʳₐ = uʳₜ + cₜₕ - uʳₕ 

    bool = xʳₐ > 1e-12 && πʳₐ > 1e-16

    return bool
end

"""
    𝝳(p::PAS, λ)

Evaluates amount of flow `δ` to shift on pas `p`.
If `δ` is less than the threshold limit of `λ` then `δ` is assumed to be zero.
"""
function 𝝳(p::PAS, λ)
    e₁, e₂, o = p.e₁, p.e₂, p.o
    
    f₁, f₂ = fₑ(e₁, o), fₑ(e₂, o)
    c₁, c₂ = cₑ(e₁), cₑ(e₂)
    c₁′, c₂′ = cₑ′(e₁), cₑ′(e₂)
    
    Δ = (c₂ - c₁)/(c₁′ + c₂′)
    
    if abs(c₂ - c₁) < λ δ = 0.0 end
    if isnan(Δ) δ = 0.0
    elseif Δ ≥ 0 δ = min(Δ, f₂)
    else δ = max(Δ, -f₁) end

    return δ
end

"""
    shift(p::PAS, δ, ϕ)

Shifts flow `δ` on pas `p`.
Argument `ϕ` determines if the arc costs need to be updated for UE or SO assignment. 
"""
function shift(p::PAS, δ)
    e₁, e₂, o = p.e₁, p.e₂, p.o
    k = o.k
    
    for a in e₁
        a.xʳ[k] += δ
        a.c = cₐ(a)
        a.c′= cₐ′(a)
    end

    for a in e₂
        a.xʳ[k] -= δ
        a.c = cₐ(a)
        a.c′= cₐ′(a)
    end
    
    return
end

"""
    MCS(a::Arc, o::Origin, G::Graph)

Develops pas for arc `a` wrt origin `o` using Maximum Cost Search method
"""
function MCS(a::Arc, o::Origin, G::Graph)
    depth, maxdepth = 1, 2
    
    i, j = a.t.v, a.h.v
    r, L₁ = o.n, o.Lᵖ
    N, A, K = G.N, G.A, G.K
    
    pᵣₕ = path(G, L₁, r, a.h)
    
    s = 1
    p = PAS(Arc[], Arc[], o)
    
    
    while depth ≤ maxdepth
        # Intialize
        l = zeros(Int64, length(N))
        for a in pᵣₕ l[a.t.v] = -1 end
        l[i] = 1
        l[j] = 1

        L₂ = Vector{Int64}(undef, length(N))
        L₂[j] = K[i,j]
        
        # Iterate
        t, h = a.t, a.h
        while true
            h = t
            
            # Maximum Cost Search
            f = 0.0
            for n in h.T 
                k = K[n,h.v]
                x = A[k].xʳ[o.k]
                c = A[k].c
                if x > 1e-12 && c > f 
                    f = c
                    t = N[n]
                    L₂[h.v] = k
                end
            end
            
            # PAS found
            if l[t.v] == -1    
                e₁ = path(G, L₁, t, a.h)
                e₂ = path(G, L₂, t, a.h)
 
                s = l[t.v]
                p = PAS(e₁, e₂, o)

                δ = 𝝳(p, 0.0)
                shift(p, δ)
                bool = ispotential(a, o, G)
                if !bool
                    return s, p
                else
                    depth += 1
                    break
                end
            # Cycle found
            elseif l[t.v] == 1
                if depth < maxdepth
                    pₕₜ = path(G, L₂, h, t)
                    if h != t push!(pₕₜ, A[K[t.v, h.v]]) end
                    δ = Inf
                    k = o.k
                    for a in pₕₜ if a.xʳ[k] ≤ δ δ = a.xʳ[k] end end
                    for a in pₕₜ 
                        a.xʳ[k] -= δ
                        a.c = cₐ(a)
                        a.c′= cₐ′(a)
                    end
                end
                depth += 1
                break
            # Continue
            else l[t.v] = 1
            end
        end
    end

    return s, p
end