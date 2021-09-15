using TrafficAssignment
using Test
using Revise
using BenchmarkTools
using DataFrames
using CSV
using Random
using Printf
Random.seed!(1104)

testnetworks = ["Anaheim"]
algorithms = [:FW, :fukushimaFW, :conjugateFW, :TAPAS]
@testset "TrafficAssignment" begin
    for network in testnetworks
        for algorithm in algorithms
            metadata, report, solution = assigntraffic(network, method=algorithm, assignment=:UE, tol=1e-10, maxiters=50, maxruntime=300, log=:off)
            refsol = DataFrame(CSV.File(joinpath(@__DIR__, "Solution\\$network.csv")))
            rg₁ = 10^report[end, :LOG₁₀RG]
            rg₂ = abs(TrafficAssignment.getrg(network, refsol))
            println("$network - $algorithm")
            @printf("   Solution relative gap : %.3e",  rg₁)
            @printf("\n   Benchmark relative gap: %.3e\n\n",  rg₂)
            @test rg₁ ≤ rg₂
        end
    end
end
# TODO: Add more test networks