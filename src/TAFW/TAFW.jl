module TAFW

using CSV
using DataFrames
using StatsBase
using Dates
using Printf

struct Node
    v::Int64                                # Node value
    T::Vector{Int64}                        # Predecessor nodes (tail node values)
    H::Vector{Int64}                        # Successor nodes (head node values)
end
# NOTE: To avoid circular reference predecessor and sucessor arrays store node indices and not the node itself

struct Origin
    n::Node                                 # Origin node
    S::Vector{Node}                         # Vector of destination nodes
    Q::Vector{Float64}                      # Vector of flows   
end

mutable struct Arc
    t::Node                                 # Tail node (From)
    h::Node                                 # Head node (To)
    V::Float64                              # Volume capacity of the arc
    d::Float64                              # Length of the arc
    tₒ::Float64                             # Free flow travel time on the arc
    α::Float64                              # BPR parameter
    β::Float64                              # BPR parameter
    x::Float64                              # Arc flow
    c::Float64                              # Arc cost
    ϕ::Bool                                 # Assignment boolean (UE => ϕ = false; SO => ϕ = true)
end

struct Graph
    name::String                            # Network name
    N::Vector{Node}                         # Vector of nodes 
    A::Vector{Arc}                          # Vector of arcs
    K::Dict{Tuple{Int64, Int64}, Int64}     # Collection of arc indices mapped to tuple of arc head and tail node value (one to one mapping - (i,j) => findfirst(x -> (x == (i,j)), [(a.t.v, a.h.v) for a in A]))
    O::Vector{Origin}                       # Vector of origins
end

include("build.jl")
include("func.jl")
include("TA.jl")

export FW, fukushimaFW, conjugateFW

end

#= ────────────────────────────────────────────────────────────────────────────────
# TODO
1. Make arc struct immutable
2. Check performance of conjugate fw againts pure fw and fukushima fw
──────────────────────────────────────────────────────────────────────────────── =#