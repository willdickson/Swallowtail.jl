
@enum StopConditionType begin
    STOP_CONDITION_TIME 
    STOP_CONDITION_STEADY
end

name_to_stop_cond = Dict(
                    "time"   => STOP_CONDITION_TIME,
                    "steady" => STOP_CONDITION_STEADY,
                   )

struct StopCondition{T}
    type::StopConditionType
    time::T
    etol::T
    function StopCondition{T}(conf::Dict{String, Any}) where {T}
        type = name_to_stop_cond[conf["type"]]
        time = conf["time"]
        etol = conf["etol"]
        return new{T}(type, time, etol)
    end
end

StopCondition(conf::Dict{String, Any}) = StopCondition{Float64}(conf)

