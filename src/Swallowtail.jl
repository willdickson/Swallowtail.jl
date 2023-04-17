module Swallowtail

using TOML
using Printf
using StaticArrays
using LinearAlgebra 
using JLD2
using Gaston


include("lattice_const.jl")
include("configuration.jl")
include("iteration.jl")
include("mesh.jl")
include("fluid.jl")
include("boundary.jl")
include("initial.jl")
include("stop.jl")
include("save.jl")
include("slbm_solver.jl")
include("viz.jl")

export runsim, plotvecs, plotvort, loaddata

function runsim(;conffile="examples/conf.toml")
    ss = SolverState(conffile)
    printconf(ss)
    runsolver(ss)
end


#if isinteractive() == false
#    runsim()
#end

end # module Swallowtail
