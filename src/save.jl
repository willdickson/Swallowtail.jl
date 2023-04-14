mutable struct SaveInfo
    count::Int
    const nstep::Int
    const npads::Int
    const directory::String
    function SaveInfo(count::Int, nstep::Int, npads::Int, directory::String)
        new(count, nstep, npads, directory)
    end
end

function SaveInfo(conf::Dict{String,Any}) 
    count = 0
    nstep = conf["nstep"]
    npads = conf["npads"]
    directory = conf["directory"]
    SaveInfo(count, nstep, npads, directory)
end

function mksave(si::SaveInfo)
    directory = savepath(si)
    mkpath(directory)
end

function savepath(si::SaveInfo)
    expanduser(si.directory)
end





