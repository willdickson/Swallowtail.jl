
struct BodyMotion{T}
end

mutable struct BodyNeighbors{T}
    num::Int
    ind::CartesianIndex
    pos::Array{T}
    u::Array{T}
end

struct Body{T}
    closed::Bool
    motion::BodyMotion{T}   
    nbrs::BodyNeighbors{T}
    pos::Array{T}
    vel::Array{T}
    ρ::Array{T}
end

struct ImmersedBoundaries{T}
end

