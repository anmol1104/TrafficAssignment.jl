"""
    itapas(G::Graph, ϵ, θ, assignment=:UE, tol=1e-5, maxIters=20, maxRunTime=600, log=:on)

improved Traffic Assignment by Paired Alternative Segments (iTAPAS) algorithm for static multi-class 
traffic assignment problem with generalized link cost function.

# Returns
a named tuple with keys `:metadata`, `:report`, and `:output`
- `metadata::String`  : Text defining the traffic assignment run 
- `report::DataFrame` : A log of total network flow, total network cost, and run time for every iteration
- `output::DataFrame` : Flow and cost for every arc from the final iteration

# Arguments
- `G::Graph`            : Network structure as `Graph`
- `assignment::Symbol`  : Assignment type; one of `:UE`, `:SO`
- `tol::Float64`        : Tolerance level for relative gap
- `maxiters::Int64`     : Maximum number of iterations
- `maxruntime::Int64`   : Maximum algorithm run time (seconds)
- `ϵ::Float64`          :
- `θ::Float64`          :
"""
function itapas(G::Graph, assignment, tol, maxiters, maxruntime, log)
    report   = DataFrame(LOG₁₀RG = Float64[], TF = Float64[], TC = Float64[], RT = Float64[])
    solution = DataFrame(FROM = Int64[], TO = Int64[], FLOW = Float64[], COST = Float64[])

    ϕ = assignment == :UE ? false : true
    
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
    for a in A a.c = cₐ(a) + ϕ * a.x * cₐ′(a) end
    for (i,o) in enumerate(O)
        r = o.n
        o.Lᵖ .= djk(G, o)
        L = o.Lᵖ
        for (j,s) in enumerate(o.S)
            qᵣₛ = o.Q[j]
            pᵣₛ = path(G, L, r, s)
            for a in pᵣₛ
                a.xʳ[i] += qᵣₛ
                a.x += qᵣₛ
                a.c = cₐ(a) + ϕ * a.x * cₐ′(a)
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
        for a in A den += a.x * a.c end
        rg = 1 - num/den
        
        # Total flow
        tf = 0.0
        for a in A tf += a.x end

        # Total cost
        tc = 0.0
        for a in A tc += a.x * a.c end

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
            L = o.Lᵖ
            T = tree(G, L)
            for (k,a) in enumerate(A)
                if a.h ∉ T && ispotential(a, o, G)
                    s, p = MCS(a, o, G)
                    if s == -1 && p ∉ P push!(P, p) end
                end
            end
            # Local shift for faster convergence
            for p in sample(P, length(P) ÷ 4) shift(p, 𝝳(p)) end
        end

        # PAS removal
        for _ in 1:40
            for (m,p) in enumerate(P)
                e₁, e₂ = p.e₁, p.e₂
                c₁, c₂ = e₁.c₁, e₂.c₂
                f₁, f₂ = e₁.f₁, e₂.f₂
                bool = (f₁ < ϵ || f₂ < ϵ) && (c₁ ≠ c₂)
                if bool deleteat!(P, m)
                else shift(p, 𝝳(p)) 
                end
            end
        end
    end
    
    for a in A
        z .= a.t.v, a.h.v, a.x, a.c
        push!(solution, z) 
    end

    metadata = "MetaData
    Network     => $(G.name)
    assignment  => $(String(assignment))
    method      => Pure Frank-Wolfe"

    return (metadata = metadata, report = report, solution = solution)
end