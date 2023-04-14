
struct FluidState{T,N}
    ρ::Array{T,N}
    u::Array{SVector{N,T},N}
    function FluidState{T,N}(dims::NTuple{N,Int}) where {T,N}
        ρ = zeros(T, dims)
        u = [@SVector zeros(T,N) for _ in CartesianIndices(dims)]
        new{T,N}(ρ,u)
    end
end

FluidState(dims::NTuple{N,Int}) where {N} = FluidState{Float64,N}(dims)
numdim(fs::FluidState{T,N}) where {T,N} = N 
Base.size(fs::FluidState) = size(fs.ρ)

function Base.zeros(fs::FluidState{T,N}) where {T,N}
    fs.ρ .= zeros(T, size(fs))
    fs.u .= [@SVector zeros(T,N) for _ in CartesianIndices(size(fs))]
    return fs
end

function zero!(fs::FluidState{T,N}) where {T,N}
    fs.ρ .= zeros(T, size(fs))
    fs.u .= [@SVector zeros(T,N) for _ in CartesianIndices(size(fs))]
end

function setstate!(fs::FluidState{T,N}, ρ::T, u::Vector{T}) where {T,N}
    length(u) == numdim(fs) || error("setstate! dimensions incompatible") 
    fs.ρ .= ρ
    fs.u .= [convert(SVector{N,T},u) for _ in CartesianIndices(size(fs))]
end

function setstate!(fs::FluidState{T,N}, ρ::T, u::SVector{N,T}) where {T,N}
    length(u) == numdim(fs) || error("setstate! dimensions incompatible") 
    fs.ρ .= ρ
    fs.u .= [u for _ in CartesianIndices(size(fs))]
end

function setstate!(fs::FluidState{T,N}, p::Array{T}, u::Array{SVector{N,T}}) where {T,N}
    fs.ρ .= ρ
    fs.u .= u
end

function setstate!(fs1::FluidState{T,N}, fs2::FluidState{T,N}) where {T,N}
    fs1.ρ .= fs2.ρ
    fs1.u .= fs2.u
end

struct Fluid{T,N}
    last::FluidState{T,N}  # Fluid state at lastious time step
    pred::FluidState{T,N}  # Fluid state at predictor step
    curr::FluidState{T,N}  # Fluid state at current time step
    function Fluid{T,N}(dims::NTuple{N,Int}) where {T,N}
        last = FluidState{T,N}(dims)
        pred = FluidState{T,N}(dims)
        curr = FluidState{T,N}(dims)
        new{T,N}(last, pred, curr)
    end
end

Fluid(dims::NTuple{N,Int}) where {N} = Fluid{Float64,N}(dims)
numdim(fs::Fluid{T,N}) where {T,N} = N 
Base.size(f::Fluid) = size(f.last.ρ)

function setstate!(f::Fluid, ρ, u)
    setstate!(f.last, ρ, u)
    setstate!(f.pred, ρ, u)
    setstate!(f.curr, ρ, u)
end

