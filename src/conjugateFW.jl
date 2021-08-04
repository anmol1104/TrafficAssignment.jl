"""
    conjugateFW(G::Graph; assignment=:UE, method=:pure, tol=1e-5, maxiters=20, maxruntime=300, log=:on)

Frank-Wolfe (conjugate) method for traffic assignment.

# Returns a named tuple with keys `:metadata`, `:report`, and `:output`
- `metadata::String`  : Text defining the traffic assignment run 
- `report::DataFrame` : A log of total network flow, total network cost, and run time for every iteration
- `output::DataFrame` : Flow and cost for every arc from the final iteration

# Arguments
- `G::Graph`                : Network structure as `Graph`
- `assignment::Symbol=:UE`  : Assignment type; one of `:UE`, `:SO`
- `tol::Float64=1e-5`       : Tolerance level for relative gap
- `maxiters::Int64=20`      : Maximum number of iterations
- `maxruntime::Int64=300`   : Maximum algorithm run time

"""
function conjugateFW(G::Graph; assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:on)
    report = DataFrame(LOG₁₀RG = Float64[], TF = Float64[], TC = Float64[])
    output = DataFrame(FROM = Int64[], TO = Int64[], FLOW = Float64[], COST = Float64[])

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
    n = 1
    for i in N for k in 1:length(A[i]) c[i][k] = cᵢⱼ(G, i, k, x[i][k], assignment) end end
    for r in R
        Lᵖ[r] = djk(G, c, r)
        for s in S[r]
            qᵣₛ = Q[r,s]
            pᵣₛ = path(Lᵖ[r], r, s)
            for (m,i) in enumerate(pᵣₛ[1:end-1])
                k = findfirst(x -> (x == pᵣₛ[m+1]), A[i])::Int64
                x[i][k] += qᵣₛ
                c[i][k] = cᵢⱼ(G, i, k, x[i][k], assignment)
            end
        end
    end
    
    # Iterate
    while true
        num , den = 0.0, 0.0
        for i in N for k in 1:length(A[i]) c[i][k] = cᵢⱼ(G, i, k, x[i][k], assignment) end end
        for r in R Lᵖ[r] = djk(G, c, r) end
        for r in R for s in S[r] num += Q[r,s] * cₑ(G, c, path(Lᵖ[r], r, s)) end end
        for i in N for k in 1:length(A[i]) den += x[i][k] * c[i][k] end end
        rg = 1 - num/den

        push!(report[!, :LOG₁₀RG], log10(abs(rg)))
        push!(report[!, :TF], sum(sum.(x)))
        push!(report[!, :TC], den)

        tₙ = now()
        runtime = (tₙ - tₒ).value/1000

        if log == :on
            if n < 10 
                @printf("\n #%02i   | %.3e | %.5e | %.5e | %.3f ", n, log10(abs(rg)), sum(sum.(x)), den, runtime)
            else 
                @printf("\n #%.0f   | %.3E | %.5E | %.5E | %.3f ", n, log10(abs(rg)), sum(sum.(x)), den, runtime) 
            end
        end

        if rg ≤ tol || n ≥ maxiters || runtime ≥ maxruntime break end

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
        λ = 0.0
        xₙ = x
        yₙ = y
        pₙ₋₁ = n == 1 ? yₙ : P[n-1]
        num, den = 0.0, 0.0
        for i in N
            for k in 1:length(A[i])
                H = derivative(x -> cᵢⱼ(G, i, k, x, assignment), x[i][k])
                num += (pₙ₋₁[i][k] - xₙ[i][k]) * H * (yₙ[i][k] - xₙ[i][k])
                den += (pₙ₋₁[i][k] - xₙ[i][k]) * H * (yₙ[i][k] - pₙ₋₁[i][k])
            end
        end
        if num/den ≤ 0.99 λ = num/den
        elseif den == 0.0 λ = 0.0
        else λ = 0.99
        end
        for i in N for k in 1:length(A[i]) p[i][k] = λ .* pₙ₋₁[i][k] + (1 - λ) .* yₙ[i][k] end end
        
        # Seach direction
        for i in N for k in 1:length(A[i]) d[i][k] = p[i][k] - x[i][k] end end
        
        # Line search
        ε = 0.01
        l = 0.0
        u = 1.0
        α = (l + u)/2
        while abs(u-l) ≥ ε
            v = 0.0
            for i in N for k in 1:length(A[i]) v += cᵢⱼ(G, i, k, x[i][k] + α * d[i][k], assignment) * d[i][k] end end
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

    for i in N
        for k in 1:length(A[i])
            push!(output[!, :FROM], i)
            push!(output[!, :TO], A[i][k])
            push!(output[!, :FLOW], x[i][k])
            push!(output[!, :COST], c[i][k])
        end
    end

    metadata =  "MetaData
    Network     => $(G.name)
    assignment  => $(String(assignment))
    method      => Conjugate Frank-Wolfe"
    
    return (metadata = metadata, report = report, output = output)
end
