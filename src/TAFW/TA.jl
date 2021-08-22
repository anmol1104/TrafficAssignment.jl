"""
    FW(G::Graph, assignment, tol, maxiters, maxruntime, log)

Frank-Wolfe method for traffic assignment.

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
"""
function FW(G::Graph, assignment, tol, maxiters, maxruntime, log)
    report   = DataFrame(LOG₁₀RG = Float64[], TF = Float64[], TC = Float64[], RT = Float64[])
    solution = DataFrame(FROM = Int64[], TO = Int64[], FLOW = Float64[], COST = Float64[])
    
    ϕ = assignment == :UE ? false : true
    
    N, A, K, O = G.N, G.A, G.K, G.O                                         # Graph  
    R  = [o.n for o in O]                                                   # Origin nodes
    y  = zeros(length(A))                                                   # Auxiliary arc flow
    p  = zeros(length(A))                                                   # Point of sight
    d  = zeros(length(A))                                                   # Search direction
    Y  = [zeros(length(A)) for _ in 1:maxiters]                             # Auxiliary link flows from each iteration (for Fukushima FW)
    P  = [zeros(length(A)) for _ in 1:maxiters]                             # Point of sight from each iteration (for Conjugate FW)
    Lᵖ = [zeros(Int64, length(N)) for _ in O]                               # Predecessor arc label (index in vector A)

    if log == :on
        print("\n iter  | LOG₁₀(RG)  | TF          | TC          |  RT (s)")
        print("\n ------|------------|-------------|-------------|--------")
    end

    # Intializate
    tₒ = now() 
    n = 0
    for a in A a.c = cₐ(a) + ϕ * a.x * cₐ′(a) end
    for (i,o) in enumerate(O)
        r = o.n
        Lᵖ[i] = djk(G, o)
        L = Lᵖ[i]
        for (j,s) in enumerate(o.S)
            qᵣₛ = o.Q[j]
            pᵣₛ = path(G, L, r, s)
            for a in pᵣₛ
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
            Lᵖ[i] = djk(G, o)
            L = Lᵖ[i]
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

        # Auxilary problem
        for (k,a) in enumerate(A) y[k] = 0.0 end
        for (i,o) in enumerate(O)
            r = o.n
            L = Lᵖ[i]
            for (j,s) in enumerate(o.S)
                qᵣₛ = o.Q[j]
                pᵣₛ = path(G, L, r, s)
                for a in pᵣₛ
                    k = K[a.t.v, a.h.v]
                    y[k] += qᵣₛ
                end
            end
        end
        
        # Point of sight
        for (k,a) in enumerate(A) p[k] = y[k] end 
        
        # Seach direction
        for (k,a) in enumerate(A) d[k] = p[k] - a.x end
        
        # Line search
        ε = 0.01
        l = 0.0
        u = 1.0
        α = (l + u)/2
        while abs(u-l) ≥ ε
            v = 0.0
            for (k,a) in enumerate(A) 
                x = a.x + α * d[k]
                v += (cₐ(a, x) + ϕ * x * cₐ′(a, x)) * d[k] 
            end
            if v < 0.0 l = α
            elseif v > 0.0 u = α
            else break
            end
            α = (l + u)/2
        end

        # Update
        for (k,a) in enumerate(A) 
            a.x += α * d[k]
            a.c = cₐ(a) + ϕ * a.x * cₐ′(a)
        end
        Y[n] = deepcopy(y)
        P[n] = deepcopy(p)
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



"""
    fukushimaFW(G::Graph, assignment, tol, maxiters, maxruntime, log)

Fukushima Frank-Wolfe method for traffic assignment.

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
"""
function fukushimaFW(G::Graph, assignment, tol, maxiters, maxruntime, log)
    report   = DataFrame(LOG₁₀RG = Float64[], TF = Float64[], TC = Float64[], RT = Float64[])
    solution = DataFrame(FROM = Int64[], TO = Int64[], FLOW = Float64[], COST = Float64[])

    ϕ = assignment == :UE ? false : true

    N, A, K, O = G.N, G.A, G.K, G.O                                         # Graph
    R  = [o.n for o in O]                                                   # Origin nodes
    y  = zeros(length(A))                                                   # Auxiliary arc flow
    p  = zeros(length(A))                                                   # Point of sight
    d  = zeros(length(A))                                                   # Search direction
    Y  = [zeros(length(A)) for _ in 1:maxiters]                             # Auxiliary link flows from each iteration (for Fukushima FW)
    P  = [zeros(length(A)) for _ in 1:maxiters]                             # Point of sight from each iteration (for Conjugate FW)
    Lᵖ = [zeros(Int64, length(N)) for _ in O]                               # Predecessor arc label (index in vector A)
    
    if log == :on
        print("\n iter  | LOG₁₀(RG)  | TF          | TC          |  RT (s)")
        print("\n ------|------------|-------------|-------------|--------")
    end

    # Intialize
    tₒ = now() 
    n = 0
    for a in A a.c = cₐ(a) + ϕ * a.x * cₐ′(a) end
    for (i,o) in enumerate(O)
        r = o.n
        Lᵖ[i] = djk(G,o)
        L = Lᵖ[i]
        for (j,s) in enumerate(o.S)
            qᵣₛ = o.Q[j]
            pᵣₛ = path(G, L, r, s)
            for a in pᵣₛ
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
            Lᵖ[i] = djk(G, o)
            L = Lᵖ[i]
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

        # Auxilary problem
        for (k,a) in enumerate(A) y[k] = 0.0 end
        for (i,o) in enumerate(O)
            r = o.n
            L = Lᵖ[i]
            for (j,s) in enumerate(o.S)
                qᵣₛ = o.Q[j]
                pᵣₛ = path(G, L, r, s)
                for a in pᵣₛ
                    k = K[a.t.v, a.h.v]
                    y[k] += qᵣₛ
                end
            end
        end

        # Point of sight
        λ = 0.5
        yₙ = y
        yₙ₋₁ = n == 1 ? yₙ : Y[n-1]
        for (k,a) in enumerate(A) p[k] = λ * yₙ₋₁[k] + (1 - λ) * yₙ[k] end

        # Seach direction
        for (k,a) in enumerate(A) d[k] = p[k] - a.x end
        
        # Line search
        ε = 0.01
        l = 0.0
        u = 1.0
        α = (l + u)/2
        while abs(u-l) ≥ ε
            v = 0.0
            for (k,a) in enumerate(A) 
                x = a.x + α * d[k]
                v += (cₐ(a, x) + ϕ * x * cₐ′(a, x)) * d[k] 
            end            
            if v < 0.0 l = α
            elseif v > 0.0 u = α
            else break
            end
            α = (l + u)/2
        end

        # Update
        for (k,a) in enumerate(A) 
            a.x += α * d[k]
            a.c = cₐ(a) + ϕ * a.x * cₐ′(a)
        end
        Y[n] = deepcopy(y)
        P[n] = deepcopy(p)
    end

    for a in A
        z .= a.t.v, a.h.v, a.x, a.c
        push!(solution, z) 
    end

    metadata =  "MetaData
    Network     => $(G.name)
    assignment  => $(String(assignment))
    method      => Fukushima FrankWolfe"
    
    return (metadata = metadata, report = report, solution = solution)
end



"""
    conjugateFW(G::Graph, assignment, tol, maxiters, maxruntime, log)

Conjugate Frank-Wolfe method for traffic assignment.

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
"""
function conjugateFW(G::Graph, assignment, tol, maxiters, maxruntime, log)
    report   = DataFrame(LOG₁₀RG = Float64[], TF = Float64[], TC = Float64[], RT = Float64[])
    solution = DataFrame(FROM = Int64[], TO = Int64[], FLOW = Float64[], COST = Float64[])

    ϕ = assignment == :UE ? false : true
    
    N, A, K, O = G.N, G.A, G.K, G.O                                         # Graph
    R  = [o.n for o in O]                                                   # Origin nodes
    y  = zeros(length(A))                                                   # Auxiliary arc flow
    p  = zeros(length(A))                                                   # Point of sight
    d  = zeros(length(A))                                                   # Search direction
    Y  = [zeros(length(A)) for _ in 1:maxiters]                             # Auxiliary link flows from each iteration (for Fukushima FW)
    P  = [zeros(length(A)) for _ in 1:maxiters]                             # Point of sight from each iteration (for Conjugate FW)
    Lᵖ = [zeros(Int64, length(N)) for _ in O]                               # Predecessor arc label (index in vector A)
    
    if log == :on
        print("\n iter  | LOG₁₀(RG)  | TF          | TC          |  RT (s)")
        print("\n ------|------------|-------------|-------------|--------")
    end

    # Intializate
    tₒ = now() 
    n = 0
    for a in A a.c = cₐ(a) + ϕ * a.x * cₐ′(a) end
    for (i,o) in enumerate(O)
        r = o.n
        Lᵖ[i] = djk(G, o)
        L = Lᵖ[i]
        for (j,s) in enumerate(o.S)
            qᵣₛ = o.Q[j]
            pᵣₛ = path(G, L, r, s)
            for a in pᵣₛ
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
            Lᵖ[i] = djk(G, o)
            L = Lᵖ[i]
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

        # Auxilary problem
        for (k,a) in enumerate(A) y[k] = 0.0 end
        for (i,o) in enumerate(O)
            r = o.n
            L = Lᵖ[i]
            for (j,s) in enumerate(o.S)
                qᵣₛ = o.Q[j]
                pᵣₛ = path(G, L, r, s)
                for a in pᵣₛ
                    k = K[a.t.v, a.h.v]
                    y[k] += qᵣₛ
                end
            end
        end

        # Point of sight
        λ = 0.0
        yₙ = y
        pₙ₋₁ = n == 1 ? yₙ : P[n-1]
        num, den = 0.0, 0.0
        for (k,a) in enumerate(A)
            num += (pₙ₋₁[k] - a.x) * cₐ′(a) * (yₙ[k] - a.x)
            den += (pₙ₋₁[k] - a.x) * cₐ′(a) * (yₙ[k] - pₙ₋₁[k])
        end
        if num/den ≤ 0.99 λ = num/den
        elseif den == 0.0 λ = 0.0
        else λ = 0.99
        end
        for (k,a) in enumerate(A) p[k] = λ .* pₙ₋₁[k] + (1 - λ) .* yₙ[k] end
        
        # Seach direction
        for (k,a) in enumerate(A) d[k] = p[k] - a.x end
        
        # Line search
        ε = 0.01
        l = 0.0
        u = 1.0
        α = (l + u)/2
        while abs(u-l) ≥ ε
            v = 0.0
            for (k,a) in enumerate(A) 
                x = a.x + α * d[k]
                v += (cₐ(a, x) + ϕ * x * cₐ′(a, x)) * d[k] 
            end
            if v < 0.0 l = α
            elseif v > 0.0 u = α
            else break
            end
            α = (l + u)/2
        end

        # Update
        for (k,a) in enumerate(A) 
            a.x += α * d[k]
            a.c = cₐ(a) + ϕ * a.x * cₐ′(a)
        end
        Y[n] = deepcopy(y)
        P[n] = deepcopy(p)
    end

    for a in A
        z .= a.t.v, a.h.v, a.x, a.c
        push!(solution, z) 
    end

    metadata =  "MetaData
    Network     => $(G.name)
    assignment  => $(String(assignment))
    method      => Conjugate Frank-Wolfe"
    
    return (metadata = metadata, report = report, solution = solution)
end


