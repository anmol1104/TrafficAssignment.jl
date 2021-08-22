using DataFrames
using CSV
using Calculus
using Plots
using Statistics

plotly()
@doc "
    traffic_assignment(;networkName, tol=0.001, maxIters=20)
"
function traffic_assignment(;networkName, tol=0.001, maxIters=20)
    networkFile = "Network\\$networkName\\network.csv"
    demandFile = "Network\\$networkName\\demand.csv"

    𝒩 = Int64[]                            # Nodes
    𝒜 = Array{Int64,1}[]                   # Arcs as adjacency list
    v₀ = Array{Float64,1}[]                 # Link volume capcity
    t₀ = Array{Float64,1}[]                 # Link free flow travel time
    α₀ = Array{Float64,1}[]                 # BPR parameters
    β₀ = Array{Float64,1}[]                 # BPR parameters

    W = Array{Int64,1}[]                    # OD pairs
    q = Float64[]                           # Demand between OD pairs
    P = Array{Array{Int64,1},1}[]           # Paths between OD pairs

    cⁱ = Float64[]                  # Stores total cost from each iteration
    xⁱₐ = Array{Array{Float64,1}}[] # Stores link flows from each iteration
    eⁱ = Float64[]                  # Stores convergent from each iteration

    # Enumerates all paths from pivot node i using nodes from a vector of open nodes 𝕏
    function enumerate(i, 𝕏)
        paths = Array{Int64,1}[]          # Vector of paths
        𝕏 = filter(x -> (x ≠ i), 𝕏)

        # Dead ends: if there is no outgoing arc from node i or if all outgoing nodes from i are closed
        if isempty(𝒜[i]) push!(paths, [i]) end
        if sum([if 𝒜[i][j] in 𝕏 1 else 0 end for j in 1:length(𝒜[i])]) == 0 push!(paths, [i]) end

        # Dveveloping paths recursively with subpaths
        for j in 1:length(𝒜[i])
            if 𝒜[i][j] in 𝕏
                for subpath in enumerate(𝒜[i][j], 𝕏)
                    push!(paths, [i])
                    append!(paths[end], subpath)
                end
            end
        end
        return paths
    end

    # Network build
    function build()
        println()
        # File import
        network = CSV.File(networkFile, types=[Int64, Int64, Float64, Float64, Float64, Float64, Float64])
        df = DataFrame(network)
        head = df[1]::Array{Int64,1}
        tail = df[2]::Array{Int64,1}
        linkcapacity = df[3]::Array{Float64,1}
        linklength = df[4]::Array{Float64,1}
        linkfft = df[5]::Array{Float64,1}
        alpha = df[6]::Array{Float64,1}
        beta = df[7]::Array{Float64,1}

        n = max(maximum(head), maximum(tail))
        for i in 1:n
            append!(𝒩, i)
            push!(𝒜, [])
            push!(v₀, [])
            push!(t₀, [])
            push!(α₀, [])
            push!(β₀, [])
        end

        for i in 1:length(head)
            append!(𝒜[head[i]], tail[i])
            append!(v₀[head[i]], linkcapacity[i])
            append!(t₀[head[i]], linkfft[i])
            append!(α₀[head[i]], alpha[i])
            append!(β₀[head[i]], beta[i])
        end

        demand = CSV.File(demandFile, types=[Int64, Int64, Float64])
        df′ = DataFrame(demand)
        origin = df′[1]::Array{Int64,1}
        destination = df′[2]::Array{Int64,1}
        demand = df′[3]::Array{Float64,1}

        for i in 1:nrow(df′)
            push!(W, [origin[i], destination[i]])
            append!(q, demand[i])
            push!(P, [])
        end

        for r in unique(origin)
            S = [W[k][2] for k in findall(x -> (x == r), origin)]
            for s in S
                index = findfirst(x -> (x == [r,s]), W)
                paths = Array{Int64,1}[]
                for path in enumerate(r, copy(𝒩))
                    if s in path
                        i = findfirst(x -> (x == s), path)
                        push!(paths, path[1:i])
                    end
                end
                append!(P[index], paths)
            end
        end
    end

    # Arc cost
    function 𝐶ₐ(i, j, x)
        α = α₀[i][j]
        β = β₀[i][j]
        t = t₀[i][j]
        v = v₀[i][j]
        return t * (1 + α * (x/v) ^ β)
    end

    # Path cost for path p given arc flows x and arc cost function cₐ
    function 𝐶ₚ(p, x, cₐ)
        cost = 0.0
        for k in 2:length(p)
            i = p[k-1]
            j = findfirst(x -> (x == p[k]), 𝒜[i])::Int64
            cost += cₐ(i, j, x[i][j])
        end
        return cost
    end

    # Djikstra's shortest path between r and s
    function sp(𝒞, r, s)
        p = Int64[]
        ℒ = [[if i==r r else -1 end, if i==r 0.0 else Inf end] for i in 𝒩]    # predecessor and cost
        𝕏 = copy(𝒩)   # Set of open nodes
        i = r
        deleteat!(𝕏, i)
        while !isempty(𝕏)
            for j in 1:length(𝒜[i])
                k = 𝒜[i][j]
                C = ℒ[i][2] + 𝒞[i][j]
                if C < ℒ[k][2] && k in 𝕏
                    ℒ[k][1] = i
                    ℒ[k][2] = C
                end
            end
            index = argmin([ℒ[i][2] for i in 𝕏])
            i = 𝕏[index]
            deleteat!(𝕏, index)
        end

        i = s
        append!(p, i)
        while i ≠ r
            i = Int(ℒ[i][1])
            append!(p, i)
        end
        reverse!(p)
        return p
    end

    # Exact Equilibration Algorithm
    function eea(cₚ)
        fₚ = Array{Float64,1}[[0.0 for _ in P[n]] for n in 1:length(W)]
        for n in 1:length(W)
            hₙ = [cₚ(n,m,0.0) for m in 1:length(P[n])]
            index = sortperm(hₙ)
            hₙ = hₙ[index]
            append!(hₙ, Inf)
            λₙ, l = 0, 0
            for k in 1:length(P[n])
                λᵏₙ = (sum(hₙ[1:k]) + q[n])/k
                if λᵏₙ > hₙ[k] && λᵏₙ ≤ hₙ[k+1]
                    l = k
                    λₙ = λᵏₙ
                    break
                end
            end
            for k′ in 1:length(P[n])
                k = findfirst(x -> (x == k′), index)
                if k ≤ l
                    fₚ[n][k′] = λₙ - hₙ[k]
                end
            end
        end
        return fₚ
    end

    # Projected Dynamical Systems: Euler method
    function pds(;cₐ=𝐶ₐ, ε=tol, kₘₐₓ=maxIters)
        fᵏₚ = [[0.0 for _ in P[n]] for n in 1:length(W)]
        xᵏₐ = Array{Float64,1}[[0.0 for j in 𝒜[i]] for i in 𝒩]

        # Intialization
        cₒ = [[cₐ(i, j, 0.0) for j in 1:length(𝒜[i])] for i in 𝒩]
        for n in 1:length(W)
            r, s = W[n]
            qᵣₛ = q[n]
            pᵣₛ = sp(cₒ, r, s)
            m = findfirst(x -> (x == pᵣₛ), P[n])
            fᵏₚ[n][m] = qᵣₛ
            for l in 2:length(pᵣₛ)
                i = pᵣₛ[l-1]
                j = findfirst(x -> (x == pᵣₛ[l]), 𝒜[i])
                xᵏₐ[i][j] += fᵏₚ[n][m]
            end
        end

        f′(i,j,x) = derivative(x -> cₐ(i,j,x))(x)
        L = sum([(f′(i, j, xᵏₐ[i][j]))^2 for i in 𝒩 for j in 1:length(𝒜[i])])^0.5
        γ = 1/L
        println("gamma: $γ")

        println("\niter: 0")
        println("total flow: ", sum(sum.(xᵏₐ)))
        println("total cost: ", sum([xᵏₐ[i][j] * cₐ(i, j, xᵏₐ[i][j]) for i in 𝒩 for j in 1:length(𝒜[i])]))
        println("convergent: NA")

        append!(cⁱ, sum([xᵏₐ[i][j] * cₐ(i, j, xᵏₐ[i][j]) for i in 𝒩 for j in 1:length(𝒜[i])]))
        push!(xⁱₐ, xᵏₐ)

        # Iteration
        cₚ = 𝐶ₚ
        converged = false
        k = 1
        while !converged
            γₖ = 1/max(k,3)
            xᵏ⁺¹ₐ = Array{Float64,1}[[0.0 for j in 𝒜[i]] for i in 𝒩]

            # Compute path flows
            c(n, m, f) = 0.5*f - (fᵏₚ[n][m] - γₖ*cₚ(P[n][m], xᵏₐ, cₐ))
            fᵏ⁺¹ₚ = eea(c)

            # Compute link flows
            for n in 1:length(W)
                for m in 1:length(P[n])
                    p = P[n][m]
                    for l in 2:length(p)
                        i = p[l-1]
                        j = findfirst(x -> (x == p[l]), 𝒜[i])
                        xᵏ⁺¹ₐ[i][j] += fᵏ⁺¹ₚ[n][m]
                    end
                end
            end

            # Convergence test
            Δ = sqrt(sum(vcat((xᵏ⁺¹ₐ - xᵏₐ)...) .^2))/sum(vcat(xᵏₐ...))
            if Δ < ε || k == kₘₐₓ converged = true end

            # Updating link flows and cost
            xᵏₐ = deepcopy(xᵏ⁺¹ₐ)

            println("\niter: $k")
            println("total flow: ", sum(sum.(xᵏₐ)))
            println("total cost: ", sum([xᵏₐ[i][j] * cₐ(i, j, xᵏₐ[i][j]) for i in 𝒩 for j in 1:length(𝒜[i])]))
            println("convergent: $Δ")

            append!(cⁱ, sum([xᵏₐ[i][j] * cₐ(i, j, xᵏₐ[i][j]) for i in 𝒩 for j in 1:length(𝒜[i])]))
            push!(xⁱₐ, xᵏₐ)
            append!(eⁱ, Δ)

            k += 1
        end
        return xᵏₐ
    end

    build()
    println(pds())
    return cⁱ, xⁱₐ, eⁱ
end
traffic_assignment(networkName="Braess-Example", maxIters=20)
