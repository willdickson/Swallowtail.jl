
#function makemesh(Δs::T, shape::Tuple{Vararg{Integer}}) where {T}
#    dims = (shape...,length(shape))
#    mesh = zeros(T,dims)
#    n = length(dims)
#    for ind in CartesianIndices(dims)
#        mesh[ind] = Δs*(ind[ind[n]]-1)
#    end
#    return mesh
#end

function meshshape(lenx::Integer, leny::Integer, Δs::T) where {T}
    nx = round(Int, lenx/Δs) + 1 
    ny = round(Int, leny/Δs) + 1
    return (nx, ny)
end

function meshshape(conf::Dict{String, Any},  Δs::T) where {T} 
    lenx = Int(conf["len"]["x"])
    leny = Int(conf["len"]["y"])
    return meshshape(lenx, leny, Δs)
end

