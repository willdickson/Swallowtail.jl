
# Schema types 
# ------------
struct NumberSchema{T<:Number}
    pos::Bool
    minval::T
    maxval::T
    required::Bool
    function NumberSchema{T}(; 
                             pos::Bool = false,
                             minval::T = typemin(T),
                             maxval::T = typemax(T),
                             required::Bool=true
                            ) where {T}
        new{T}(pos, minval, maxval, required)
    end
end


struct DirectorySchema
    exist::Bool
    required::Bool
    function DirectorySchema(;exist::Bool=false, required::Bool=true)
        new(exist, required)
    end
end


struct StringOptSchema
    options::Array{String}
    required::Bool
    function StringOptSchema(;options::Array{String}, required::Bool=true)
        new(options, required)
    end
end

# Schema creator
# --------------
function makeschema(T::DataType=Float64)

    # mesh
    lenschema = Dict(
                     "x" => NumberSchema{T}(pos=true), 
                     "y" => NumberSchema{T}(pos=true),
                    )
    meshschema = Dict(
                      "ds"  => NumberSchema{T}(pos=true), 
                      "len" => lenschema,
                     )

    # save
    saveschema = Dict(
                      "nstep" => NumberSchema{Int}(pos=true),
                      "npads" => NumberSchema{Int}(pos=true),
                      "directory" => DirectorySchema(),
                     )
    # fluid
    fluidschema = Dict(
                       "kvisc" => NumberSchema{T}(pos=true),
                       "density" => NumberSchema{T}(pos=true),
                      )
    veloschema = Dict(
                      "x" => NumberSchema{T}(),
                      "y" => NumberSchema{T}(),
                     )
    # init
    initopt = ["constant",]
    initschema = Dict(
                      "type" => StringOptSchema(options=initopt),
                      "velocity" => veloschema,
                     )
    # bndry
    bndryopt = ["inflow", "outflow", "moving", "noslip", "slip"]
    sideschema = Dict( 
                      "type" => StringOptSchema(options=bndryopt), 
                      "value" => NumberSchema{T}(),
                     )

    bndryschema = Dict(
                       "left" => sideschema,
                       "right" => sideschema, 
                       "top" => sideschema,
                       "bottom" => sideschema,
                      )


    # stop
    stopopt = ["time", "steady"]
    stopschema = Dict(
                      "type" => StringOptSchema(options=stopopt),
                      "time" => NumberSchema{T}(pos=true),
                      "etol" => NumberSchema{T}(pos=true),
                     )

    # conf schema
    confschema = Dict(
        "mesh"  => meshschema, 
        "save"  => saveschema, 
        "fluid" => fluidschema, 
        "init"  => initschema, 
        "bndry" => bndryschema, 
        "stop"  => stopschema
       )
    return confschema
end


# Load conf from toml file
function loadconf(filename::String, T::DataType=Float64)
    conf = Dict()
    try
        conf = TOML.parsefile(expanduser(filename))
    catch e
        error("unable to parse file: $filename, $e")
    end
    return conf
end


# Check conf to ensure conforms with schema
function checkconf(conf::Dict, T::DataType=Float64)
    schema = makeschema(T)
    ok, confmod, msg = checkconf("root", conf, schema)
    if !ok
        msg = "conf$msg" 
    end
    return ok, confmod, msg
end


function checkconf(k::String,  conf::Dict, schema::Dict)
    ok, msg = true, ""
    confmod = Dict{String, Any}()
    for (k,v) in schema 
        if k in keys(conf)
            ok, val, msg = checkconf(k, conf[k], v)
            if !ok
                msg = "[\"$k\"]$msg"
                break
            end
            confmod[k] = val
        else
            if v["required"]
                ok = false
                msg = " missing key $k"
                break
            end
        end
    end
    return ok, confmod, msg 
end


function checkconf(k::String, val, schema::NumberSchema{T}) where {T}
    ok, msg = true, ""
    valmod = convert(T, val)
    if valmod > schema.maxval
        ok = false
        msg = "[\"$k\"] must be ≤ maxval"
        return ok, valmod, msg
    end
    if valmod < schema.minval
        ok = false
        msg = "[\"$k\"] must be ≥ minval"
        return ok, valmod, msg
    end
    if schema.pos && valmod ≤ 0
        ok = false
        msg = "[\"$k\"] must be > 0"
        return ok, valmod, msg
    end
    return ok, valmod, msg
end


function checkconf(k::String, val, schema::StringOptSchema)
    ok, msg = true, ""
    valmod = string(val)
    if valmod ∉ schema.options
        ok = false
        msg = "[\"$k\"] not in options"
    end
    return ok, valmod, msg
end


function checkconf(k::String, val, schema::DirectorySchema)
    ok, msg = true, ""
    valmod = string(val)
    if schema.exist && not isdir(valmod)
        ok = false
        msg = "[\"$k\"] directory must exist"
    end
    return ok, val, msg
end

