
mutable struct Iteration{T}
    t::T     # Current time in simulation 
    n::Int   # Current interation number
    function Iteration{T}(t::T=0.0,n::Int=0) where {T}
        new{T}(t,n)
    end
end

Interation(t=0.0, n=0) = Iteration{Float64}(t,n)

function next!(it::Iteration{T}, Δt::T) where {T}
    it.t += Δt
    it.n += 1
end

