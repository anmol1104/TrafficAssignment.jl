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

    ğ© = Int64[]                            # Nodes
    ğ = Array{Int64,1}[]                   # Arcs as adjacency list
    vâ = Array{Float64,1}[]                 # Link volume capcity
    tâ = Array{Float64,1}[]                 # Link free flow travel time
    Î±â = Array{Float64,1}[]                 # BPR parameters
    Î²â = Array{Float64,1}[]                 # BPR parameters

    W = Array{Int64,1}[]                    # OD pairs
    q = Float64[]                           # Demand between OD pairs
    P = Array{Array{Int64,1},1}[]           # Paths between OD pairs

    câ± = Float64[]                  # Stores total cost from each iteration
    xâ±â = Array{Array{Float64,1}}[] # Stores link flows from each iteration
    eâ± = Float64[]                  # Stores convergent from each iteration

    # Enumerates all paths from pivot node i using nodes from a vector of open nodes ğ
    function enumerate(i, ğ)
        paths = Array{Int64,1}[]          # Vector of paths
        ğ = filter(x -> (x â  i), ğ)

        # Dead ends: if there is no outgoing arc from node i or if all outgoing nodes from i are closed
        if isempty(ğ[i]) push!(paths, [i]) end
        if sum([if ğ[i][j] in ğ 1 else 0 end for j in 1:length(ğ[i])]) == 0 push!(paths, [i]) end

        # Dveveloping paths recursively with subpaths
        for j in 1:length(ğ[i])
            if ğ[i][j] in ğ
                for subpath in enumerate(ğ[i][j], ğ)
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
            append!(ğ©, i)
            push!(ğ, [])
            push!(vâ, [])
            push!(tâ, [])
            push!(Î±â, [])
            push!(Î²â, [])
        end

        for i in 1:length(head)
            append!(ğ[head[i]], tail[i])
            append!(vâ[head[i]], linkcapacity[i])
            append!(tâ[head[i]], linkfft[i])
            append!(Î±â[head[i]], alpha[i])
            append!(Î²â[head[i]], beta[i])
        end

        demand = CSV.File(demandFile, types=[Int64, Int64, Float64])
        dfâ² = DataFrame(demand)
        origin = dfâ²[1]::Array{Int64,1}
        destination = dfâ²[2]::Array{Int64,1}
        demand = dfâ²[3]::Array{Float64,1}

        for i in 1:nrow(dfâ²)
            push!(W, [origin[i], destination[i]])
            append!(q, demand[i])
            push!(P, [])
        end

        for r in unique(origin)
            S = [W[k][2] for k in findall(x -> (x == r), origin)]
            for s in S
                index = findfirst(x -> (x == [r,s]), W)
                paths = Array{Int64,1}[]
                for path in enumerate(r, copy(ğ©))
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
    function ğ¶â(i, j, x)
        Î± = Î±â[i][j]
        Î² = Î²â[i][j]
        t = tâ[i][j]
        v = vâ[i][j]
        return t * (1 + Î± * (x/v) ^ Î²)
    end

    # Path cost for path p given arc flows x and arc cost function câ
    function ğ¶â(p, x, câ)
        cost = 0.0
        for k in 2:length(p)
            i = p[k-1]
            j = findfirst(x -> (x == p[k]), ğ[i])::Int64
            cost += câ(i, j, x[i][j])
        end
        return cost
    end

    # Djikstra's shortest path between r and s
    function sp(ğ, r, s)
        p = Int64[]
        â = [[if i==r r else -1 end, if i==r 0.0 else Inf end] for i in ğ©]    # predecessor and cost
        ğ = copy(ğ©)   # Set of open nodes
        i = r
        deleteat!(ğ, i)
        while !isempty(ğ)
            for j in 1:length(ğ[i])
                k = ğ[i][j]
                C = â[i][2] + ğ[i][j]
                if C < â[k][2] && k in ğ
                    â[k][1] = i
                    â[k][2] = C
                end
            end
            index = argmin([â[i][2] for i in ğ])
            i = ğ[index]
            deleteat!(ğ, index)
        end

        i = s
        append!(p, i)
        while i â  r
            i = Int(â[i][1])
            append!(p, i)
        end
        reverse!(p)
        return p
    end

    # General Equilibration Algorithm
    function gea(;câ, Îµ, kâââ)
        câ = ğ¶â
        Îµâ² = 0.000001
        kâ²âââ = 20

        xáµâ = Array{Float64,1}[[0.0 for j in ğ[i]] for i in ğ©]
        fáµâ = Array{Float64,1}[[0.0 for _ in P[n]] for n in 1:length(W)]

        # Intialization
        cáµâ = [[câ(i, j, 0.0) for j in 1:length(ğ[i])] for i in ğ©]
        for n in 1:length(W)
            r, s = W[n]
            qáµ£â = q[n]
            páµ£â = sp(cáµâ, r, s)
            m = findfirst(x -> (x == páµ£â), P[n])
            fáµâ[n][m] = qáµ£â
            for l in 2:length(páµ£â)
                i = páµ£â[l-1]
                j = findfirst(x -> (x == páµ£â[l]), ğ[i])
                xáµâ[i][j] += fáµâ[n][m]
            end
        end

        # Iteration
        allConverged = false
        k = 1
        while !allConverged
            xáµâºÂ¹â = deepcopy(xáµâ)
            fáµâºÂ¹â = deepcopy(fáµâ)

            for n in 1:length(W)
                odConverged=false
                kâ² = 1
                while !odConverged
                    c = [câ(P[n][m], xáµâºÂ¹â, câ) for m in 1:length(P[n])]
                    indices = findall(x -> (x > 0), fáµâºÂ¹â[n])

                    # Active path with highest costs
                    mâ = indices[argmax([c[m] for m in indices])]
                    # Path with least cost
                    mâ = argmin(c)

                    X = [[P[n][mâ][l-1], P[n][mâ][l]] for l in 2:length(P[n][mâ])]
                    Y = [[P[n][mâ][l-1], P[n][mâ][l]] for l in 2:length(P[n][mâ])]
                    Z = intersect(X, Y)
                    num = c[mâ] - c[mâ]
                    den = Float64(length(X) + length(Y) - 2*length(Z))
                    if num == den == 0 Î´ = 0.0
                    else Î´ = min(num/den, fáµâºÂ¹â[n][mâ])
                    end

                    # OD convergence test
                    if Î´ < Îµâ² || kâ² == kâ²âââ odConverged = true end

                    # Updation
                    fáµâºÂ¹â[n][mâ] -= Î´
                    fáµâºÂ¹â[n][mâ] += Î´
                    for l in 2:length(P[n][mâ])
                        i = P[n][mâ][l-1]
                        j = findfirst(x -> (x == P[n][mâ][l]), ğ[i])
                        xáµâºÂ¹â[i][j] -= Î´
                    end
                    for l in 2:length(P[n][mâ])
                        i = P[n][mâ][l-1]
                        j = findfirst(x -> (x == P[n][mâ][l]), ğ[i])
                        xáµâºÂ¹â[i][j] += Î´
                    end

                    kâ² += 1
                end
            end

            # Overall convergence test
            num = 0.0
            den = 0.0
            for i in ğ©
                for j in 1:length(ğ[i])
                    num += (xáµâºÂ¹â[i][j] - xáµâ[i][j])^2
                    den += xáµâ[i][j]
                end
            end
            Î = sqrt(num)/den
            if Î < Îµ || k == kâââ allConverged = true end

            # Updation
            xáµâ = deepcopy(xáµâºÂ¹â)
            fáµâ = deepcopy(fáµâºÂ¹â)

            k += 1
        end
        return xáµâ
    end

    # Variational Inequality: Projection method
    function vi(;câ=ğ¶â, Îµ=tol, kâââ=maxIters)
        xáµâ = Array{Float64,1}[[0.0 for j in ğ[i]] for i in ğ©]

        # Intialization
        câ = [[câ(i, j, 0.0) for j in 1:length(ğ[i])] for i in ğ©]
        for n in 1:length(W)
            r, s = W[n]
            qáµ£â = q[n]
            páµ£â = sp(câ, r, s)
            for l in 2:length(páµ£â)
                i = páµ£â[l-1]
                j = findfirst(x -> (x == páµ£â[l]), ğ[i])
                xáµâ[i][j] += qáµ£â
            end
        end

        fâ²(i,j,x) = derivative(x -> câ(i,j,x))(x)
        L = sum([(fâ²(i, j, xáµâ[i][j]))^2 for i in ğ© for j in 1:length(ğ[i])])^0.5
        Î³ = 1/L
        println("gamma: $Î³")

        println("\niter: 0")
        println("total flow: ", sum(sum.(xáµâ)))
        println("total cost: ", sum([xáµâ[i][j] * câ(i, j, xáµâ[i][j]) for i in ğ© for j in 1:length(ğ[i])]))
        println("convergent: NA")

        append!(câ±, sum([xáµâ[i][j] * câ(i, j, xáµâ[i][j]) for i in ğ© for j in 1:length(ğ[i])]))
        push!(xâ±â, xáµâ)

        # Iteration
        converged = false
        k = 1
        while !converged
            xáµâºÂ¹â = Array{Float64,1}[[0.0 for j in ğ[i]] for i in ğ©]

            # Auxilary problem
            c(i,j,x) = x - (xáµâ[i][j] - Î³*câ(i,j,xáµâ[i][j]))
            xáµâºÂ¹â = gea(câ=c, Îµ=0.000001, kâââ=20)

            # Convergence test
            Î = sqrt(sum(vcat((xáµâºÂ¹â - xáµâ)...) .^2))/sum(vcat(xáµâ...))
            if Î < Îµ || k == kâââ converged = true end

            # Updating link flows and cost
            xáµâ = deepcopy(xáµâºÂ¹â)

            println("\niter: $k")
            println("total flow: ", sum(sum.(xáµâ)))
            println("total cost: ", sum([xáµâ[i][j] * câ(i, j, xáµâ[i][j]) for i in ğ© for j in 1:length(ğ[i])]))
            println("convergent: $Î")

            append!(câ±, sum([xáµâ[i][j] * câ(i, j, xáµâ[i][j]) for i in ğ© for j in 1:length(ğ[i])]))
            push!(xâ±â, xáµâ)
            append!(eâ±, Î)

            k += 1
        end
        return xáµâ
    end

    build()
    println(vi())
    return câ±, xâ±â, eâ±
end
traffic_assignment(networkName="Braess-Example", tol=0.001, maxIters=20)
