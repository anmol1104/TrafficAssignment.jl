# Auxilary functions
# Arc functions ────────────────────────────────────────────────────────────────────────────────
"""
    cₐ(i, k, x, method)

Returns arc cost for arc `(i,j)` on graph `G` given arc flow `x` (`j = A[i][k]`)
"""
function cₐ(G::Graph, i, k, x, assignment)
    α = G.α[i][k]
    β = G.β[i][k]
    tₒ= G.T[i][k]
    V = G.V[i][k]
    
    t = tₒ * (1 + α * (abs(x)/V) ^ β)

    if assignment == :UE
        c = t
    end
    if assignment == :SO
        t′= β == 0 ? 0.0 : tₒ * α * β * (abs(x) ^ (β - 1))/(V ^ β)
        c = t + x * t′
    end
    return c
end

"""
    c′ₐ(i, k, x, method)

Returns derivative of arc cost for arc `(i,j)` on graph `G` given arc flow `x` (`j = A[i][k]`)
"""
function cₐ′(G::Graph, i, k, x, assignment)
    α = G.α[i][k]
    β = G.β[i][k]
    tₒ= G.T[i][k]
    V = G.V[i][k]

    t = tₒ * (1 + α * (abs(x)/V) ^ β)
    t′= β == 0 ? 0.0 : tₒ * α * β * (abs(x) ^ (β - 1))/(V ^ β)

    if assignment == :UE
        c′ = t′
    end
    if assignment == :SO
        t″ = β == 0 || β == 1 ?  0.0 : tₒ * α * β * (β - 1) * (abs(x) ^ (β - 2))/(V ^ β)
        c′ = 2t′+ x * t″ 
    end
    return c′
end



# Segment functions ────────────────────────────────────────────────────────────────────────────────
"""
    cₑ(G::Graph, c, e)

Returns cost for segment `e` on graph `G` given arc costs `c`
"""
function cₑ(G::Graph, c, e)
    A = G.A 
    v = 0.0
    for (n,i) in enumerate(e[1:end-1])
        j = e[n+1]
        k = findfirst(x -> (x == j), A[i])::Int64
        v += c[i][k]
    end
    return v
end

"""
    fₑ(G::Graph, x, e)

Returns minimum flow on segment `e` on graph `G`  given arc flows `x`
"""
function fₑ(G::Graph, x, e)
    A = G.A
    f = zeros(length(e)-1)
    for (n,i) in enumerate(e[1:end-1])
        j = e[n+1]
        k = findfirst(x -> (x == j), A[i])::Int64
        f[n] = x[i][k]
    end
    return minimum(f)
end


# Djikstra functions ────────────────────────────────────────────────────────────────────────────────
"""
    djk(G::Graph, c, r)

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


# ────────────────────────────────────────────────────────────────────────────────
# Checks if arc a fails reduced cost optimal conditions for origin r
function ispotential(a, r)
    i, j = a
    k = findfirst(x -> (x == j), A[i])::Int64
    pᵣᵢ = path(Lᵖ[r], r, i)
    pᵣⱼ = path(Lᵖ[r], r, j)
    uʳᵢ = cₑ(G, c, pᵣᵢ)
    uʳⱼ = cₑ(G, c, pᵣⱼ)
    cʳᵢⱼ = cₐ[m][i][k]
    πʳᵢⱼ = uʳᵢ + cʳᵢⱼ - uʳⱼ
    if xʳₐ[r][i][k] > ϵ && πʳᵢⱼ > θ return (true, πʳᵢⱼ)
    else return (false, 0.0) end
end

# Checks if PAS p assosciated with origin rₒ can be eliminated
function isbad(p, rₒ)
    e₁, e₂ = p
    c₁, c₂ = cₑ(G, c, e₁), cₑ(G, c, e₂)
    f₁, f₂ = fₑ(G, xʳ[rₒ], e₁), fₑ(G, xʳ[rₒ], e₂)
    if (f₁ < ϵ || f₂ < ϵ) && (c₁ ≠ c₂) return true
    else return false end
end

# Shifts flows from higher cost segment to lower cost segment of PAS p 
# on its assosciated origin rₒ, given cost difference is greater than λ
function shift(p, rₒ, λ)
    e₁, e₂ = p
    m = Mᵣ[rₒ]

    c₁, c₂ = cₑ(e₁, cₐ[m]), cₑ(e₂, cₐ[m])
    if abs(c₂ - c₁) < λ return end

    c′₁, c′₂ = cₑ(e₁, c′ₐ[m]), cₑ(e₂, c′ₐ[m])
    f₁, f₂ = fₑ(e₁, xʳₐ[rₒ]), fₑ(e₂, xʳₐ[rₒ])
    Δ = (c₂ - c₁)/(c′₁ + c′₂)
    if isnan(Δ) δ = 0.0
    elseif Δ ≥ 0 δ = min(Δ, f₂)
    else δ = max(Δ, -f₁) end

    for (n,i) in enumerate(e₁[1:end-1])
        j = e₁[n+1]
        k = findfirst(x -> (x == j), A[i])::Int64
        xʳₐ[rₒ][i][k] += δ
        xₐ[i][k] += δ
        for m in M
            cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k], assignment)
            c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k], assignment)
        end
    end
    for (n,i) in enumerate(e₂[1:end-1])
        j = e₂[n+1]
        k = findfirst(x -> (x == j), A[i])::Int64
        xʳₐ[rₒ][i][k] -= δ
        xₐ[i][k] -= δ
        for m in M
            cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k], assignment)
            c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k], assignment)
        end
    end
end

# PAS identification for arc a wrt origin r using Maximum Cost Search
function MCS(a, r)
    depth, maxdepth = 1, 2
    flag = true

    i, j = a
    e₁, e₂ = Int64[], Int64[]
    pᵣⱼ = path(Lᵣ[r], r, j)

    while flag
        # Intialize
        lₖ = [if k ∈ a 1 elseif k ∉ pᵣⱼ 0 else -1 end for k in N]
        L = [if k == j i else -1 end for k in N]

        # Iterate
        t = i
        while true
            h = t

            f = 0.0
            for p in A′[h]
                k = findfirst(x -> (x == h), A[p])
                x = xʳₐ[r][p][k]
                c = cᵢⱼ(p, k, m, x, assignment)
                if x > ϵ && c > f f, t = c, p end
            end

            L[h] = t
            if lₖ[t] == -1      # PAS found
                e₁ = path(Lᵣ[r], t, j)
                e₂ = path(L, t, j)
                shift((e₁, e₂), r, 0)
                bool,_ = ispotential(a, r)
                if !bool || depth == maxdepth flag = false
                else depth += 1 end
                break
            elseif lₖ[t] == 1   # Cycle found
                if depth == maxdepth flag = false
                else
                    if h == t pₕₜ = Int64[]
                    else
                        pₕₜ = path(L, h, t)
                        push!(pₕₜ, h)
                        δ = fₑ(pₕₜ, xʳₐ[r])
                    end
                    for (n,i) in enumerate(pₕₜ[1:end-1])
                        k = findfirst(x -> (x == pₕₜ[n+1]), A[i])::Int64
                        xʳₐ[r][i][k] -= δ
                        xₐ[i][k] -= δ
                        for m in M
                            cₐ[m][i][k] = cᵢⱼ(i, k, m, xₐ[i][k], assignment)
                            c′ₐ[m][i][k] = c′ᵢⱼ(i, k, m, xₐ[i][k], assignment)
                        end
                    end
                    depth += 1
                end
                break
            else                # Continue
                lₖ[t] = 1
            end
        end
    end
    p = (e₁, e₂)
    return p
end