
struct LatticeConst{T}
    cs  ::T                         # Lattice speed of sound
    a1  ::T                         # 1st equilib function constant 1/(cs^2)
    a2  ::T                         # 2nd equilib function constant 1/(2*cs^4)
    a3  ::T                         # 3rd equilib function constant 1/(2*cs^2)
    w   ::Array{T}                  # Array of lattice weights
    e   ::Vector{SVector{2, Int}}   # Array of lattice vectors
    ind ::Array{CartesianIndex{2}}  # Index offsets associated with lattice vectors
    num ::Int                       # Number of lattice vectors
    function LatticeConst{T}() where {T}
        cs = 1/sqrt(3)
        a1 = 1/(cs^2)
        a2 = 1/(2*cs^4)
        a3 = 1/(2*cs^2)
        w  = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36] 
        e  = [ 
              SA[ 0, 0], 
              SA[ 1, 0], 
              SA[ 0, 1], 
              SA[-1, 0], 
              SA[ 0,-1], 
              SA[ 1, 1], 
              SA[-1, 1], 
              SA[-1,-1], 
              SA[ 1,-1], 
             ] 
        ind  = [CartesianIndex(v...) for v in e]
        num  = 9
        new{T}(cs, a1, a2, a3, w, e, ind, num)
    end
end

LatticeConst() = LatticeConst{Float64}()



