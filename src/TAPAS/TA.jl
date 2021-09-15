"""
    itapas(G::Graph, tol, maxIters, maxRunTime, log)

improved Traffic Assignment by Paired Alternative Segments (iTAPAS) algorithm for static multi-class 
traffic assignment problem with generalized link cost function.

# Returns
a named tuple with keys `:metadata`, `:report`, and `:output`
- `metadata::String`  : Text defining the traffic assignment run 
- `report::DataFrame` : A log of total network flow, total network cost, and run time for every iteration
- `output::DataFrame` : Flow and cost for every arc from the final iteration

# Arguments
- `G::Graph`            : Network structure as `Graph`
- `tol::Float64`        : Tolerance level for relative gap
- `maxiters::Int64`     : Maximum number of iterations
- `maxruntime::Int64`   : Maximum algorithm run time (seconds)
"""
function itapas(G::Graph, tol, maxiters, maxruntime, log)
    report   = DataFrame(LOG₁₀RG = Float64[], TF = Float64[], TC = Float64[], RT = Float64[])
    solution = DataFrame(FROM = Int64[], TO = Int64[], FLOW = Float64[], COST = Float64[])
    
    N, A, K, O = G.N, G.A, G.K, G.O                                         # Graph
    R  = [o.n for o in O]                                                   # Origin nodes
    P  = PAS[]                                                              # PASs
    
    if log == :on
        print("\n iter  | logRG      | TF          | TC          | RT (s) ")
        print("\n ------|------------|-------------|-------------|--------")
    end
    
    # Intialize
    tₒ = now() 
    n = 0
    for a in A 
        a.c = cₐ(a)
        a.c′= cₐ′(a)
    end 
    for (i,o) in enumerate(O)
        r = o.n
        o.Lᵖ .= djk(G, o)
        L = o.Lᵖ
        for (j,s) in enumerate(o.S)
            qᵣₛ = o.Q[j]
            pᵣₛ = path(G, L, r, s)
            for a in pᵣₛ
                a.xʳ[i] += qᵣₛ
                a.c = cₐ(a)
            end
        end
    end

    # Iterate
    z = zeros(4)
    while true
        # Relative gap
        num, den = 0.0, 0.0
        for (i,o) in enumerate(O)
            r = o.n
            o.Lᵖ .= djk(G, o)
            L = o.Lᵖ
            for (j,s) in enumerate(o.S)
                qᵣₛ = o.Q[j]
                pᵣₛ = path(G, L, r, s)
                for a in pᵣₛ num += qᵣₛ * a.c end
            end
        end
        for a in A den += sum(a.xʳ) * a.c end
        rg = 1 - num/den
        
        # Total flow
        tf = 0.0
        for a in A tf += sum(a.xʳ) end

        # Total cost
        tc = 0.0
        for a in A tc += sum(a.xʳ) * a.c end

        # Run time
        tₙ = now()
        runtime = (tₙ - tₒ).value/1000
        
        z .= log10(abs(rg)), tf, tc, runtime
        push!(report, z)
        
        if log == :on @printf("\n #%02i   | %.3e | %.5e | %.5e | %.3f ", n, z...) end

        if rg ≤ tol || n ≥ maxiters || runtime ≥ maxruntime break end

        n += 1

        # Indentify potential arc => Find PAS for this arc => Perform flow shift on this PAS
        for (i,o) in enumerate(O)
            o.Lᵖ .= djk(G, o)
            L = o.Lᵖ
            T = tree(G, L)
            for (k,a) in enumerate(A)
                if a.h ∈ T[a.t.v] continue end
                bool = ispotential(a, o, G)
                if bool
                    s, p = MCS(a, o, G)
                    if s == -1 && p ∉ P push!(P, p) end
                end
            end
            # Local shift for faster convergence
            for p in sample(P, length(P) ÷ 4) 
                δ = 𝝳(p, rg/1000)
                shift(p, δ)  
            end
        end
        
        # PAS removal
        for _ in 1:40
            for (m,p) in enumerate(P)
                e₁, e₂, o = p.e₁, p.e₂, p.o
                f₁, f₂ = fₑ(e₁, o), fₑ(e₂, o)
                c₁, c₂ = cₑ(e₁), cₑ(e₂)
                bool = (f₁ < 1e-12 || f₂ < 1e-12) && (c₁ ≠ c₂)
                if bool deleteat!(P, m)
                else shift(p, 𝝳(p, rg/1000)) 
                end
            end
        end
    end
    
    for a in A
        z .= a.t.v, a.h.v, sum(a.xʳ), a.c
        push!(solution, z) 
    end

    assignment = A[begin].ϕ == false ? :UE : :SO
    
    metadata = "MetaData
    Network     => $(G.name)
    assignment  => $(String(assignment))
    method      => Pure Frank-Wolfe"

    return (metadata = metadata, report = report, solution = solution)
end