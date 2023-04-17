
function vorticity2d(u, Δs)
    ind = CartesianIndices(u)

    ux = [u[i][1] for i in ind]
    uy = [u[i][2] for i in ind]

    ix = CartesianIndex((1,0))
    iy = CartesianIndex((0,1))
    
    ∂x(f, i, Δx) = 0.5*(f[i+ix] - f[i-ix])/Δx
    ∂y(f, i, Δy) = 0.5*(f[i+iy] - f[i-iy])/Δy

    ω = zeros(size(ux))
    n, m = size(ux)

    for i in CartesianIndices((2:n-1,2:m-1))
        ω[i] = ∂x(uy,i,Δs) - ∂y(ux,i,Δs)
    end
    return ω
end


function loaddata(datadir::String, savenum::Int)
    npad = lpad(savenum,6,"0")
    datadir = expanduser(datadir)
    filename = "$npad.jld2"
    filepath = joinpath(datadir, filename)
    println("loading $filepath")
    return load_object(filepath)
end


function plotvort(datadir::String, savenum::Int)
    ss = loaddata(datadir, savenum)
    curr = ss.fluid.curr
    ind = CartesianIndices(size(curr.ρ))
    Δs = ss.Δs
    scale = 1.0
    ω = vorticity2d(curr.u, Δs)
    n, m = size(ω)
    x = [Δs*i for i in 1:n]
    y = [Δs*i for i in 1:m]
    heatmap(x, y, ω', Axes(cbrange=(-0.3,0.3), size="square"))
end

function plotvecs(datadir::String, savenum::Int)

    ss = loaddata(datadir, savenum)
    curr = ss.fluid.curr
    Δs = ss.Δs
    ind = CartesianIndices(size(curr.ρ))
    scale = 1.0
    x = vec([Δs*i[1] for i in ind])
    y = vec([Δs*i[2] for i in ind])
    ux = scale*vec([curr.u[i][1] for i in ind])
    uy = scale*vec([curr.u[i][2] for i in ind])

    k = 1
    plot(
         x[1:k:end], 
         y[1:k:end], 
         supp=[ux[1:k:end] uy[1:k:end]], 
         w=:vectors, 
         Axes(size="square"),
        )

end

function plotvecs(datadir::String, savenum::Int, sectnum::Int)
    ss = loaddata(datadir, savenum)
    curr = ss.fluid.curr
    Δs = ss.Δs
    ind = CartesianIndices(size(curr.ρ))
    scale = 1.0
    x = vec([Δs*i[1] for i in ind][sectnum,:])
    y = vec([Δs*i[2] for i in ind][sectnum,:])
    ux = scale*vec([curr.u[i][1] for i in ind][sectnum,:])
    uy = scale*vec([curr.u[i][2] for i in ind][sectnum,:])

    k = 1
    plot(
         x[1:k:end], 
         y[1:k:end], 
         supp=[ux[1:k:end] uy[1:k:end]], 
         w=:vectors, 
         Axes(size="square"),
        )

end

function plotvecs(datadir::String, savenum::Int, sectnums::Array{Int})
    ss = loaddata(datadir, savenum)
    curr = ss.fluid.curr
    Δs = ss.Δs
    ind = CartesianIndices(size(curr.ρ))
    scale = 1.0
    x = vec([Δs*i[1] for i in ind][sectnums,:])
    y = vec([Δs*i[2] for i in ind][sectnums,:])
    ux = scale*vec([curr.u[i][1] for i in ind][sectnums,:])
    uy = scale*vec([curr.u[i][2] for i in ind][sectnums,:])

    k = 1
    plot(
         x[1:k:end]*0.0, 
         y[1:k:end], 
         supp=[ux[1:k:end] uy[1:k:end]], 
         w=:vectors, 
         Axes(size="square"),
        )

end

