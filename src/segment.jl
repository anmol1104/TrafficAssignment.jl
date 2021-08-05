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

"""
    cₑ(G, e)

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
    fₑ(G, e)

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