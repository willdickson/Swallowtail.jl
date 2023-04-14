
@enum InitialConditionType begin 
    INIT_CONST 
    INIT_VECT_X 
    INIT_VECT_Y 
    INIT_FIELD
end

struct InitialCondition{T} 
    type::InitialConditionType
    u::Array{T}
end

function InitialCondition{T}(
                             conf::Dict{String, Any},          
                             shape::Tuple{Vararg{Integer}},
                            ) where {T}
    conf_str = conf["type"]
    if conf_str == "constant"
        type = INIT_CONST
        u = zeros(T,2)
        u[1] = convert(T, conf["velocity"]["x"])
        u[2] = convert(T, conf["velocity"]["y"])
        InitialCondition{T}(type, u)
    else
        error("intial cond type=$conf_str not implemeted yet")
    end
end

function InitialCondition(conf::Dict{String, Any}, shape::Tuple{Vararg{Integer}})
    InitialCondition{Float64}(conf, shape)
end

