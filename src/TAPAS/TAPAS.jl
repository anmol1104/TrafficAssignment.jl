module TrafficAssignment

using CSV
using DataFrames
using StatsBase
using Dates
using Printf
using Plots

struct Node
    v::Int64                                # Node value
    T::Vector{Int64}                        # Predecessor nodes (tail node values)
    H::Vector{Int64}                        # Successor nodes (head node values)
end
# NOTE: To avoid circular reference predecessor and sucessor arrays store node indices and not the node itself

struct Origin
    n::Node                                 # Origin node
    k::Int64                                # Origin node index in vector G.O; k = findfirst(x -> (x == o), O)
    S::Vector{Node}                         # Vector of destination nodes
    Q::Vector{Float64}                      # Vector of flows   
    Lᵖ::Vector{Int64}                       # Predecessor label
end

mutable struct Arc
    t::Node                                 # Tail node (From)
    h::Node                                 # Head node (To)
    V::Float64                              # Volume capacity of the arc
    d::Float64                              # Length of the arc
    tₒ::Float64                             # Free flow travel time on the arc
    α::Float64                              # BPR parameter
    β::Float64                              # BPR parameter
    xʳ::Vector{Float64}                     # Origin based arc flow
    c::Float64                              # Arc cost
end

mutable struct Segment
    s::Vector{Arc}                          # Segment
    f::Vector{Float64}                      # Minimum origin-based flow on the segment
    c::Float64                              # Sum of arc cost on the segment
    c′::Float64                             # Sum of derivative of arc cost on the segment
end

struct PAS
    e₁::Segment                             # First segment
    e₂::Segment                             # Second segment
    o::Origin                               # Assosciated origin
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

export itapas

end

#= ────────────────────────────────────────────────────────────────────────────────
# TODO
1. Make arc and segment struct immutable
2. Test tapas against benchmarks from Xie, Nie and Liu (2018) - A greedy path based algorithm
──────────────────────────────────────────────────────────────────────────────── =#