"""
    FW(G::Graph, assignment, tol, maxiters, maxruntime, log)

Frank-Wolfe method for traffic assignment.

# Returns a named tuple with keys `:metadata`, `:report`, and `:output`
- `metadata::String`  : Text defining the traffic assignment run 
- `report::DataFrame` : A log of total network flow, total network cost, and run time for every iteration
- `output::DataFrame` : Flow and cost for every arc from the final iteration

# Arguments
- `G::Graph`            : Network structure as `Graph`
- `assignment::Symbol`  : Assignment type; one of `:UE`, `:SO`
- `tol::Float64`        : Tolerance level for relative gap
- `maxiters::Int64`     : Maximum number of iterations
- `maxruntime::Int64`   : Maximum algorithm run time
"""
function FW(G::Graph, assignment, tol, maxiters, maxruntime, log)
    report   = DataFrame(LOG₁₀RG = Float64[], TF = Float64[], TC = Float64[], RT = Float64[])
    solution = DataFrame(FROM = Int64[], TO = Int64[], FLOW = Float64[], COST = Float64[])
    
    N, A, R, S, Q = G.N, G.A, G.R, G.S, G.Q                                         # Graph  
    x = [zeros(length(A[i])) for i in N]                                            # Arc flow
    c = [zeros(length(A[i])) for i in N]                                            # Arc cost
    y = [zeros(length(A[i])) for i in N]                                            # Auxiliary arc flow
    p = [zeros(length(A[i])) for i in N]                                            # Point of sight
    d = [zeros(length(A[i])) for i in N]                                            # Search direction
    Y = Array{Array{Float64,1},1}[]                                                 # Stores auxiliary link flows from each iteration (for Fukushima FW)
    P = Array{Array{Float64,1},1}[]                                                 # Stores point of sight from each iteration (for Conjugate FW)
    Lᵖ= Dict(r => [if i==r r else -1 end for i in N] for r in R)                    # Djikstra's oredecessor labels for every origin node
    
    if log == :on
        print("\n iter  | LOG₁₀(RG)  | TF          | TC          |  RT (s)")
        print("\n ------|------------|-------------|-------------|--------")
    end

    # Intialization
    tₒ = now() 
    n = 0
    for i in N for k in 1:length(A[i]) c[i][k] = cₐ(G, i, k, x[i][k], assignment) end end
    for r in R
        Lᵖ[r] = djk(G, c, r)
        for s in S[r]
            qᵣₛ = Q[r,s]
            pᵣₛ = path(Lᵖ[r], r, s)
            for (m,i) in enumerate(pᵣₛ[1:end-1])
                k = findfirst(x -> (x == pᵣₛ[m+1]), A[i])::Int64
                x[i][k] += qᵣₛ
                c[i][k] = cₐ(G, i, k, x[i][k], assignment)
            end
        end
    end
    
    # Iterate
    while true
        num , den = 0.0, 0.0
        for i in N for k in 1:length(A[i]) c[i][k] = cₐ(G, i, k, x[i][k], assignment) end end
        for r in R Lᵖ[r] = djk(G, c, r) end
        for r in R for s in S[r] num += Q[r,s] * cₑ(G, c, path(Lᵖ[r], r, s)) end end
        for i in N for k in 1:length(A[i]) den += x[i][k] * c[i][k] end end
        rg = 1 - num/den
        
        tₙ = now()
        runtime = (tₙ - tₒ).value/1000
        
        push!(report, [log10(abs(rg)), sum(sum.(x)), den, runtime])
        
        if log == :on
            @printf("\n #%02i   | %.3e | %.5e | %.5e | %.3f ", n, log10(abs(rg)), sum(sum.(x)), den, runtime)
        end

        if rg ≤ tol || n + 1 ≥ maxiters || runtime ≥ maxruntime break end

        n += 1

        # Auxilary problem
        for i in N for k in 1:length(A[i]) y[i][k] = 0.0 end end
        for r in R
            for s in S[r]
                qᵣₛ = Q[r,s]
                pᵣₛ = path(Lᵖ[r], r, s)
                for (m,i) in enumerate(pᵣₛ[1:end-1])
                    k = findfirst(x -> (x == pᵣₛ[m+1]), A[i])::Int64
                    y[i][k] += qᵣₛ
                end
            end
        end

        # Point of sight
        for i in N for k in 1:length(A[i]) p[i][k] = y[i][k] end end
        
        # Seach direction
        for i in N for k in 1:length(A[i]) d[i][k] = p[i][k] - x[i][k] end end
        
        # Line search
        ε = 0.01
        l = 0.0
        u = 1.0
        α = (l + u)/2
        while abs(u-l) ≥ ε
            v = 0.0
            for i in N for k in 1:length(A[i]) v += cₐ(G, i, k, x[i][k] + α * d[i][k], assignment) * d[i][k] end end
            if v < 0.0 l = α
            elseif v > 0.0 u = α
            else break
            end
            α = (l + u)/2
        end
        
        # Update
        for i in N for k in 1:length(A[i]) x[i][k] += α * d[i][k] end end
        push!(Y, y)
        push!(P, p)
    end

    for i in N for k in 1:length(A[i]) push!(solution, [i, A[i][k], x[i][k], c[i][k]]) end end

    metadata = "MetaData
    Network     => $(G.name)
    assignment  => $(String(assignment))
    method      => Pure Frank-Wolfe"
    
    return (metadata = metadata, report = report, solution = solution)
end
