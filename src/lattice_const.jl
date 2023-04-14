
struct LatticeConst{T}
    cs  ::T                         # Speed of sound
    a1  ::T                         # Equilib function constant 1/cs
    a2  ::T                         # Equilib function constant 1/(2*cs^4)
    a3  ::T                         # Equilib function constant 1/(2*cs^2)
    w   ::Array{T}                  # Array of lattice weights
    e   ::Vector{SVector{2, Int}}   # Array of lattice vectors
    ind ::Array{CartesianIndex{2}}  # Lattice vectors index offsets
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



