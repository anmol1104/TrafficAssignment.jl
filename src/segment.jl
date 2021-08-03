"""
    cₑ(G, e)

Returns cost for segment `e` on graph `G`
"""
function cₑ(G, e)
    _, A, C, = G
    c = 0.0
    for (n,i) in enumerate(e[1:end-1])
        j = e[n+1]
        k = findfirst(x -> (x == j), A[i])::Int64
        c += C[i][k]
    end
    return c
end

"""
    fₑ(G, e)

Returns minimum flow on segment `e` on graph `G`
"""
function fₑ(G, e)
    _, A, _, X = G
    f = zeros(length(e)-1)
    for (n,i) in enumerate(e[1:end-1])
        j = e[n+1]
        k = findfirst(x -> (x == j), A[i])::Int64
        f[n] = X[i][k]
    end
    return minimum(f)
end