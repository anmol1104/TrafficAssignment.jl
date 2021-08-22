module TrafficAssignment

include("TAFW\\TAFW.jl")
include("TAPAS\\TAPAS.jl")
include("main.jl")

using .TAFW

export assigntraffic

end
#= ────────────────────────────────────────────────────────────────────────────────
# TODO:
1. Finish jobs in FW
2. Complete TAPAS
3. Complete tests
──────────────────────────────────────────────────────────────────────────────── =#