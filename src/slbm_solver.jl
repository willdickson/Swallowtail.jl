
struct SolverState{T}
    conf::Dict{String,Any}      # Configuration dictionary
    iter::Iteration{T}          # Current iteration information
    lbm::LatticeConst{T}        # Lattice Boltzam constants
    Δt::T                       # Time step
    Δs::T                       # Mesh spacing
    ρ0::T                       # Reference density
    ν ::T                       # Kinematic viscosity
    τ ::T                       # Relaxation parameter
    fluid::Fluid{T,2}           # Fluid states (last, pred, curr)
    bndry::DomainBoundaries{T}  # Domain boundary conditions
    initc::InitialCondition{T}  # Simulation initial condition
    stopc::StopCondition{T}     # Simulation stop condition
    savei::SaveInfo             # Save infomation
    
    function SolverState{T}(conf::Dict{String, Any}) where {T}
        ok, conf, msg = checkconf(conf)
        ok || error(msg)
        iter = Iteration{T}()
        lbm = LatticeConst{T}()
        Δt = conf["mesh"]["ds"]
        Δs = conf["mesh"]["ds"]
        ρ0 = conf["fluid"]["density"]
        ν  = conf["fluid"]["kvisc"]
        τ  = T(0.5) + ν/(Δt*lbm.cs^2)
        shape = meshshape(conf["mesh"], Δs) 
        fluid = Fluid{T,2}(shape)
        bndry = DomainBoundaries{T}(conf["bndry"], shape)
        initc = InitialCondition{T}(conf["init"], shape)
        stopc = StopCondition{T}(conf["stop"])
        savei = SaveInfo(conf["save"])
        new{T}(conf, iter, lbm, Δt, Δs, ρ0, ν, τ, fluid, bndry, initc, stopc, savei) 
    end
end


Base.size(ss::SolverState) = size(ss.fluid)
numdim(ss::SolverState) = numdim(ss.fluid)
SolverState(conf::Dict{String, Any}) = SolverState{Float64}(conf)


function SolverState(conffile::String) 
    conf = loadconf(conffile)
    SolverState(conf)
end


function runsolver(ss::SolverState{T}) where {T}
    println()
    println("running slbm solver")
    println("------------------------")
    println()
    initialize!(ss)     # Set initial condition
    setboundary!(ss)    # Set boundary condition
    done::Bool = false  # Flag indicating completion 
    while !done
        next!(ss.iter, ss.Δt)
        predictor!(ss)
        corrector!(ss)
        ibcorrector!(ss)
        savedata(ss)
        done = isdone(ss)
        #printinfo(ss)
    end
end


function initialize!(ss::SolverState)
    mksave(ss.savei)
    if ss.initc.type == INIT_CONST
        setstate!(ss.fluid, ss.ρ0, ss.initc.u)
    else
        error("initial cond type = $(ss.initc.type) not implemented")
    end
end


function setboundary!(ss::SolverState{T}) where {T}
    setboundary!(ss.fluid.last, ss.bndry)
    setboundary!(ss.fluid.pred, ss.bndry)
    setboundary!(ss.fluid.curr, ss.bndry)
end

function setboundary!(fs::FluidState{T,N}, bndry::DomainBoundaries{T}) where {T,N}
    setboundary!(fs.ρ, fs.u, bndry)
end


function setboundary!(
        ρ::Array{T,N}, 
        u::Array{SVector{N,T}}, 
        bndry::DomainBoundaries{T},
    ) where {T,N}
    for b in [bndry.left, bndry.right, bndry.top, bndry.bottom]
        if b.type in [BOUNDARY_INFLOW, BOUNDARY_MOVING, BOUNDARY_NOSLIP]
            u[b.inds] .= [b.u for _ in b.inds]
        elseif b.type == BOUNDARY_OUTFLOW
            u[b.inds] .= u[b.adj1]
        elseif b.type == BOUNDARY_SLIP
            u[b.inds] .= [dot(u[ind],b.para)*b.para for ind in b.adj1]
        else
            error("setboundary!: unknown boundary type $(b.type)")
        end
        ρ[b.inds] .= (T(4.0).*ρ[b.adj1] .- ρ[b.adj2])/T(3.0)
    end
end


function predictor!(ss::SolverState{T}) where {T}
    predictor!( 
               fluid  = ss.fluid,
               evals  = ss.lbm.e,  
               wvals  = ss.lbm.w,  
               einds  = ss.lbm.ind, 
               a1     = ss.lbm.a1, 
               a2     = ss.lbm.a2, 
               a3     = ss.lbm.a3,
               ) 
    setboundary!(ss.fluid.pred, ss.bndry)
end

function predictor!(;
        fluid::Fluid{T,N},
        evals::Vector{SVector{N,Int}},
        wvals::Vector{T},
        einds::Vector{CartesianIndex{N}},
        a1::T,
        a2::T,
        a3::T,
    ) where {T,N}

    ρ_last = fluid.last.ρ
    u_last = fluid.last.u
    ρ_pred = fluid.pred.ρ
    u_pred = fluid.pred.u
    ρ_curr = fluid.curr.ρ
    u_curr = fluid.curr.u

    n, m = size(ρ_last)
    interior_inds = CartesianIndices((2:n-1,2:m-1)) # NOTE: 2D specific
    lattice_vals = zip(evals, wvals, einds)

    Threads.@threads for i in interior_inds
    #for i in interior_inds 
        ρ_last[i] = ρ_curr[i]
        u_last[i] = u_curr[i]
        ρ_pred[i] = 0.0
        u_pred[i] = @SVector zeros(T,N)
        for (ek, wk, ik) in lattice_vals 
            feq = equilib(
                          ρ  = ρ_curr[i-ik], 
                          u  = u_curr[i-ik], 
                          w  = wk, 
                          e  = ek, 
                          a1 = a1, 
                          a2 = a2, 
                          a3 = a3,
                         )
            ρ_pred[i] += feq
            u_pred[i] += feq*ek
        end
        u_pred[i] = u_pred[i]/ρ_pred[i]
    end

end


function corrector!(ss::SolverState{T}) where {T}
    corrector!( 
               fluid  = ss.fluid,
               evals  = ss.lbm.e,  
               wvals  = ss.lbm.w,  
               einds  = ss.lbm.ind, 
               a1     = ss.lbm.a1, 
               a2     = ss.lbm.a2, 
               a3     = ss.lbm.a3,
               τ      = ss.τ,
               ) 
    setboundary!(ss.fluid.curr, ss.bndry)
end


function corrector!(;
        fluid::Fluid{T,N},
        evals::Vector{SVector{N,Int}},
        wvals::Vector{T},
        einds::Vector{CartesianIndex{N}},
        a1::T,
        a2::T,
        a3::T,
        τ ::T,
    ) where {T,N}

    ρ_last = fluid.last.ρ
    u_last = fluid.last.u
    ρ_pred = fluid.pred.ρ
    u_pred = fluid.pred.u
    ρ_curr = fluid.curr.ρ
    u_curr = fluid.curr.u

    n, m = size(ρ_last)
    interior_inds = CartesianIndices((2:n-1,2:m-1)) # NOTE: 2D specific
    lattice_vals = zip(evals, wvals, einds)

    Threads.@threads for i in interior_inds
    #for i in interior_inds 
        u_curr[i] = ρ_pred[i]*u_pred[i]
        for (ek, wk, ik) in lattice_vals
            feq = equilib(
                          ρ  = ρ_pred[i+ik], 
                          u  = u_pred[i+ik], 
                          w  = wk, 
                          e  = ek, 
                          a1 = a1, 
                          a2 = a2, 
                          a3 = a3,
                         )
            u_curr[i] += (τ - T(1.0))*feq*ek
        end
        ρ_curr[i]  = ρ_pred[i]
        u_curr[i] -= (τ - T(1.0))*ρ_last[i]*u_last[i]
        u_curr[i] /= ρ_curr[i]
    end
end


@inline function equilib(;
        ρ::T,              # fluid density
        u::SVector{2,T},   # fluid velocity
        w::T,              # lattice vector weight
        e::SVector{2,Int}, # lattice vector
        a1::T,             # equilib function constant 1
        a2::T,             # equilib function constant 2
        a3::T,             # equilib function constant 3
    ) where {T}
    uu = sum(u.*u)
    eu = sum(e.*u)
    eu2 = eu*eu
    return ρ*w*(T(1.0) + a1*eu + a2*eu2 - a3*uu)
end


function ibcorrector!(ss::SolverState{T}) where{T}

end


function savedata(ss::SolverState) 
    if mod(ss.iter.n, ss.savei.nstep) == 0
        ss.savei.count += 1
        countstr = lpad(ss.savei.count, ss.savei.npads, "0")
        filename = "$countstr.jld2" 
        @printf(
                "saving: %s, count = %d, n = %d, t =  %1.3f\n", 
                filename,
                ss.savei.count, 
                ss.iter.n,
                ss.iter.t, 
               )
        filepath = joinpath(savepath(ss.savei), filename)
        save_object(filepath, ss)
    end
end


function printinfo(ss::SolverState) 
    @printf "n = %s, t = %1.3f\n" lpad(ss.iter.n, 8 , "0") ss.iter.t 
end


function isdone(ss::SolverState) 
    done = false
    if ss.stopc.type == STOP_CONDITION_TIME
        if ss.iter.t >= ss.stopc.time
            done = true
        end
    elseif ss.stopc.type == STOP_CONDITION_STEADY 
        done = isconverged(ss,ss.stopc.etol)
    end
    return done
    #return true 
end


function isconverged(ss::SolverState, etol::T) where {T}
    done = false
    nx, ny = size(ss)
    d = numdim(ss)
    ind = CartesianIndices((2:nx-1, 2:ny-1, 1:d))
    mag_last = map(norm, ss.last.u[ind])
    mag_curr = map(norm, ss.curr.u[ind])
    max_rerr = maximum(abs.(mag_last .- mag_curr)./mag_last)
    if max_rerr < etol
        done = true
    end
    return done
end


function printconf(ss::SolverState{T}) where {T}
    println()
    println("simulation configuration")
    println("------------------------")
    println("T    = ", T)
    println("dim  = ", numdim(ss))
    println("size = ", size(ss))
    println("Δt   = ", ss.Δt)
    println("Δs   = ", ss.Δs)
    println("ρ0   = ", ss.ρ0)
    println("ν    = ", ss.ν)
    bndrydict = Dict(
                     "bndry.left"   => ss.bndry.left, 
                     "bndry.right"  => ss.bndry.right,
                     "bndry.top"    => ss.bndry.top,
                     "bndry.bottom" => ss.bndry.bottom,
                    )
    for (name, bndry) in bndrydict
        println(name)
        println(" type = ", bndry.type)
        if bndry.type in [BOUNDARY_INFLOW, BOUNDARY_MOVING]
            println(" value = ", bndry.u[1,1])
        end
    end
    println("initc")
    if size(ss.initc.u) == (2,)
        println(" type = ", ss.initc.type)
        println(" u    = ", ss.initc.u)
    else
        println(" type    = ", ss.initc.type)
        println(" size(u) = ", size(ss.initc.u))
    end
    println("stopc")
    println(" type = ", ss.stopc.type)
    println(" time = ", ss.stopc.time)
    if !isnan(ss.stopc.etol)
        println(" etol = ", ss.stopc.etol)
    end
    println(" niter = ", floor(Int,ss.stopc.time/ss.Δt))
    println("save")
    println(" nstep = ", ss.savei.nstep)
    println(" npads = ", ss.savei.npads)
    println(" directory = ", ss.savei.directory)
end

