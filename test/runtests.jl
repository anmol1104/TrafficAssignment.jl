using TrafficAssignment
using Test
using Revise
using CSV
using DataFrames
using BenchmarkTools

assigntraffic("Anaheim", method=:FW, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:on)
@benchmark assigntraffic("Anaheim", method=:FW, assignment=:UE, tol=1e-5, maxiters=20, maxruntime=300, log=:off)

#TrafficAssignment.compare("Anaheim", methods=[:FW, :fukushimaFW, :conjugateFW])

#=
testnetworks = ["Anaheim", "SiouxFalls"]
@testset "FW.jl" begin
    for network in testnetworks
        metadata, report, solution = assigntraffic(network, method=:FW, assignment=:UE, tol=1e-10, maxiters=50, maxruntime=300, log=:off)
        refsol = DataFrame(CSV.File(joinpath(@__DIR__, "Solution\\$network.csv")))
        rg₁ = 10^report[end, :LOG₁₀RG]
        rg₂ = abs(TrafficAssignment.getrg(network, refsol))
        @test rg₁ ≤ rg₂
    end
end
=#