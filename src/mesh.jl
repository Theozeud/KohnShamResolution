#=
    Mesh structure contains sorted points

=#

struct Mesh{T} 

    points::Vector{T}   

    function Mesh(points::AbstractVector{T})  where T <: Real
        _points = union(sort(points))
        new{T}(_points)
    end
end

@inline Base.eltype(::Mesh{T}) where T = T
@inline Base.eachindex(m::Mesh) = eachindex(m.points)
@inline Base.firstindex(m::Mesh) = firstindex(m.points)
@inline Base.lastindex(m::Mesh) = lastindex(m.points)
@inline Base.getindex(m::Mesh, n::Int) = m.points[n]
@inline Base.getindex(m::Mesh, ur::UnitRange{Int64}) = m.points[ur]
@inline Base.setindex!(m::Mesh{T}, val::T, n::Int) where T = m.points[n] = val
@inline Base.first(m::Mesh) = m[firstindex(m)]
@inline Base.last(m::Mesh) = m[lastindex(m)]
@inline Base.length(m::Mesh) = length(m.points)
@inline Base.size(m::Mesh) = size(m.points)

@inline function findindex(m, x)
    if x ≤ m[end]
        return searchsortedlast(m.points, x)
    else
        return lastindex(m)+1
    end
end

   
@inline Base.iterate(m::Mesh, state = 1) = state > length(m) ? nothing : (m[state],state+1)

######################################
# LINEAR MESH
######################################

linmesh(a, b, n; T = Float64) = Mesh(T.(LinRange(a,b,n)))

######################################
# GEOMETRIC MESH
######################################

function geometricrange(a,b,n; T = Float64, s)
    R = zeros(T,n)
    R[1] = a
    hn = (one(T)-T(s))/(one(T) - T(s)^n)*(b-a)
    H = zeros(T,n-1)
    H[end] = hn
    for i ∈ n-1:-1:2
        H[i-1]= T(s) * H[i]
    end 
    for i ∈ 2:n
        R[i] = R[i-1] + H[i-1]
    end
    R
end

geometricmesh(a,b,n; T = Float64, s) = Mesh(geometricrange(a,b,n; T = T, s = s))