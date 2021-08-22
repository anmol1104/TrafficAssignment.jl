# Auxilary functions
# Arc functions ────────────────────────────────────────────────────────────────────────────────
"""
    cₐ(a::Arc, x=a.x)

Returns arc cost for `a::Arc` for arc flow `x`
"""
function cₐ(a::Arc, x=a.x)
    tₒ= a.tₒ
    α = a.α
    β = a.β
    V = a.V

    t = tₒ * (1 + α * (x/V) ^ β)

    c = t
    return c
end

"""
    cₐ′(a::Arc, x=a.x)

Returns first derivative of arc cost wrt arc flow for `a::Arc` at arc flow `x`
"""
function cₐ′(a::Arc, x=a.x)
    tₒ= a.tₒ
    α = a.α
    β = a.β
    V = a.V

    t′= β == 0 ? 0.0 : tₒ * α * β * (x ^ (β - 1))/(V ^ β)

    c′ = t′
    return c′
end


"""
    cₐ″(a::Arc, x=a.x)

Returns second derivative of arc cost wrt arc flow for `a::Arc` at arc flow `x`
"""
function cₐ″(a::Arc, x=a.x)
    tₒ= a.tₒ
    α = a.α
    β = a.β
    x = a.x
    V = a.V

    t″ = β == 0 || β == 1 ?  0.0 : tₒ * α * β * (β - 1) * (x ^ (β - 2))/(V ^ β)

    c″ = t″
    return c″
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
path(G::Graph, o::Origin, s::Node) = path(G, o.Lᵖ, o.n, s)  # Shorter version

# ────────────────────────────────────────────────────────────────────────────────


# Identfies if arc a is a potential arc wrt flow from origin o, where k = findfirst(x -> (x == o), O)
function ispotential(a::Arc, o::Origin, G::Graph)
    k = o.k
    xʳₐ = a.xʳ[k]
    cₜₕ = a.c

    t = a.t
    pᵣₜ = path(G, o, t)
    uʳₕ = 0.0
    for a in pᵣₜ uʳₕ += a.c end
    
    h = a.h
    pᵣₕ = path(G, o, h)
    uʳₜ = 0.0
    for a in pᵣₕ uʳₜ += a.c end
    
    πʳₐ = uʳₜ + cₜₕ - uʳₕ 

    bool = xʳₐ > 1e-12 && πʳₐ > 1e-16

    return bool
end

# Shifts flows from higher cost segment to lower cost segment of PAS p 
# on its assosciated origin rₒ, given cost difference is greater than λ
function 𝝳(p::PAS)#, λ)
    e₁, e₂, o = p.e₁, p.e₂, p.o
    k = o.k

    f₁ = e₁.f[k]
    c₁ = e₁.c
    c₁′= e₁.c′

    f₂ = e₂.f[k]
    c₂ = e₂.c
    c₂′= e₂.c′

    #if abs(c₂ - c₁) < λ return 0.0 end
    
    Δ = (c₂ - c₁)/(c₁′ + c₂′)
    if isnan(Δ) δ = 0.0
    elseif Δ ≥ 0 δ = min(Δ, f₂)
    else δ = max(Δ, -f₁) end
    
    return δ
end

function shift(p::PAS, δ)
    # TODO: Adjust cost vals for assingnment

    e₁, e₂, o = p.e₁, p.e₂, p.o
    k = o.k
    
    f, c, c′ = Inf, 0.0, 0.0
    for a in e₁.s
        a.xʳ[k] += δ
        a.x += δ
        a.c = cₐ(a)

        if a.xʳ[k] < f f = a.xʳ[k] end
        c += a.c
        c′+= cₐ′(a)
    end
    e₁.f[k], e₁.c, e₁.c′ = f, c, c′
    
    f, c, c′ = Inf, 0.0, 0.0
    for a in e₂.s
        a.xʳ[k] += δ
        a.x += δ
        a.c = cₐ(a)

        if a.xʳ[k] < f f = a.xʳ[k] end
        c += a.c
        c′+= cₐ′(a)
    end
    e₂.f[k], e₂.c, e₂.c′ = f, c, c′
    
    return
end


# PAS identification for arc a wrt origin r using Maximum Cost Search
function MCS(a::Arc, o::Origin, G::Graph)
    depth, maxdepth = 1, 2
    flag = true

    N, A = G.N, G.A

    L₁ = o.Lᵖ
    t, h = a.t, a.h
    i, j = t.v, h.v
    e₁, e₂ = Segment(Arc[], Inf, 0.0, 0.0), Segment(Arc[], Inf, 0.0, 0.0)
    p = PAS(e₁, e₂, o)

    pᵣⱼ = path(G, L₁, r, j)

    while flag
        # Intialize
        l = zeros(Int64, length(N))
        for a in pᵣⱼ l[a.t] = -1 end
        l[i] = 1
        l[j] = 1

        L₂ = Vector{Node}(undef, length(N))
        L₂[j] = t 


        # Iterate
        v = i
        while true
            
            f = 0.0
            for u in N[v].T 
                k = findfirst(x -> (x == N[v].v), N[u].H)::Int64
                x = xʳₐ[r][u][k]
                k = K[u,v]
                a = A[k]
                c = a.c
                if x > 1e-12 && c > f f, t = c, p end
            end
            L₂[v] = N[u]
            v = u   
            
            # PAS found
            if l[v] == -1      
                s₁ = path(G, L₁, t, j)
                s₂ = path(G, L₂, t, j)
                shift(p, 𝝳(p))
                bool,_ = ispotential(a, o, G)
                if !bool || depth == maxdepth flag = false
                else depth += 1 end
                break
            # Cycle found
            elseif l[v] == 1
                if depth == maxdepth flag = false
                else
                    pₕₜ = path(G, L₂, h, t)
                    push!(pₕₜ, a)   
                    δ = Inf
                    k = o.k
                    for a in pₕₜ if a.xʳ[k] ≤ δ δ = a.xʳ[k] end end
                    shift(PAS(pₕₜ , Int64[], o), δ)
                    depth += 1
                end
                break
            # Continue
            else l[v] = 1
            end
        end
    end

    return p
end