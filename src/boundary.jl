
@enum BoundaryType begin
    BOUNDARY_INFLOW 
    BOUNDARY_OUTFLOW 
    BOUNDARY_MOVING 
    BOUNDARY_NOSLIP 
    BOUNDARY_SLIP
end

name_to_bndry_cond = Dict( 
    "inflow"  => BOUNDARY_INFLOW, 
    "outflow" => BOUNDARY_OUTFLOW, 
    "moving"  => BOUNDARY_MOVING, 
    "noslip"  => BOUNDARY_NOSLIP, 
    "slip"    => BOUNDARY_SLIP,
   )

struct BoundaryCondition{T}
    type::BoundaryType
    u::Vector{T}
    inds::CartesianIndices
    adj1::CartesianIndices
    adj2::CartesianIndices
    para::Vector{T}
    norm::Vector{T}
    function BoundaryCondition{T}(
            type::BoundaryType, 
            u::Vector{T},
            inds::CartesianIndices, 
            adj1::CartesianIndices,
            adj2::CartesianIndices,
            para::Vector{T},
            norm::Vector{T},
        )  where {T}
        new{T}(type, u, inds, adj1, adj2, para, norm) 
    end
end

function BoundaryCondition(
        type::BoundaryType, 
        u::Vector, 
        inds::CartesianIndices, 
        adj1::CartesianIndices,
        adj2::CartesianIndices,
        para::Vector,
        norm::Vector,
    )  
    BoundaryCondition{Float64}(type, u, inds, adj1, adj2, para, norm) 
end

function bcond_left(T::DataType, conf::Dict{String, Any}, shape::Tuple{Vararg{Integer}})
    para = T[0, 1]
    norm = T[1, 0]
    type = name_to_bndry_cond[conf["type"]] 
    value = boundary_value(conf)
    u = T[value, 0]
    nx, ny = shape
    if type == BOUNDARY_INFLOW
        yrange = 2:ny-1
        u = value*norm
    elseif type == BOUNDARY_MOVING
        yrange = 1:ny
        u = value*para
    else
        yrange = 1:ny
        u =T[0.0, 0.0]
    end
    inds = CartesianIndices((1:1, yrange))
    adj1 = CartesianIndices((2:2, yrange))
    adj2 = CartesianIndices((3:3, yrange))
    BoundaryCondition{T}(type, u, inds, adj1, adj2, para, norm)
end

function bcond_right(T::DataType, conf::Dict{String, Any}, shape::Tuple{Vararg{Integer}}) 
    para = T[0,1]
    norm = T[1,0]
    type = name_to_bndry_cond[conf["type"]] 
    value = boundary_value(conf)
    u = T[value, 0]
    nx, ny = shape
    if type == BOUNDARY_INFLOW
        yrange = 2:ny-1
        u = value*norm
    elseif type == BOUNDARY_MOVING
        yrange = 1:ny
        u = value*para
    else
        yrange = 1:ny
        u =T[0.0, 0.0]
    end
    inds = CartesianIndices((nx-0:nx-0, yrange))
    adj1 = CartesianIndices((nx-1:nx-1, yrange))
    adj2 = CartesianIndices((nx-2:nx-2, yrange))
    BoundaryCondition{T}(type, u, inds, adj1, adj2, para, norm)
end

function bcond_top(T::DataType, conf::Dict{String, Any}, shape::Tuple{Vararg{Integer}})
    para = T[1,0]
    norm = T[0,1]
    type = name_to_bndry_cond[conf["type"]] 
    value = boundary_value(conf)
    nx, ny = shape
    if type == BOUNDARY_INFLOW
        xrange = 2:nx-1
        u = value*norm
    elseif type == BOUNDARY_MOVING
        xrange = 1:nx
        u = value*para
    else
        xrange = 1:nx
        u =T[0.0, 0.0]
    end
    inds = CartesianIndices((xrange, ny-0:ny-0))
    adj1 = CartesianIndices((xrange, ny-1:ny-1))
    adj2 = CartesianIndices((xrange, ny-2:ny-2))
    BoundaryCondition{T}(type, u, inds, adj1, adj2, para, norm)
end

function bcond_bottom(T::DataType, conf::Dict{String, Any}, shape::Tuple{Vararg{Integer}})
    para = T[1,0]
    norm = T[0,1]
    type = name_to_bndry_cond[conf["type"]] 
    value = boundary_value(conf)
    u = T[0, value]
    nx, ny = shape
    if type == BOUNDARY_INFLOW
        xrange = 2:nx-1
        u = value*norm
    elseif type == BOUNDARY_MOVING
        xrange = 1:nx
        u = value*para
    else
        xrange = 1:nx
        u =T[0.0, 0.0]
    end
    inds = CartesianIndices((xrange, 1:1))
    adj1 = CartesianIndices((xrange, 2:2))
    adj2 = CartesianIndices((xrange, 3:3))
    BoundaryCondition{T}(type, u, inds, adj1, adj2, para, norm)
end

struct DomainBoundaries{T}
    left::BoundaryCondition{T}
    right::BoundaryCondition{T}
    top::BoundaryCondition{T}
    bottom::BoundaryCondition{T}
    function DomainBoundaries{T}(conf::Dict{String, Any}, shape::Tuple{Vararg{Integer}},) where {T}
        left = bcond_left(T, conf["left"], shape)
        right = bcond_right(T, conf["right"], shape)
        top = bcond_top(T, conf["top"], shape)
        bottom = bcond_bottom(T, conf["bottom"], shape)
        new{T}(left, right, top, bottom)
    end
end

DomainBoundaries(conf::Dict{String, Any}) = DomainBoundaries{Float64}(conf)


function boundary_value(conf::Dict{String, Any})
    value = conf["value"]
    if conf["type"] == BOUNDARY_NOSLIP
        value = zero(value)
    end
    return value
end

