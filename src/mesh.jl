abstract type AbstractMesh end

struct OneDMesh{T} <: AbstractMesh
    points::Vector{T}
end

function mesh(points::AbstractArray{T})  where T<: Real
    _points = union(sort(points))
    OneDMesh(_points)
end

linmesh(a, b, n, T = Float64) = mesh(T.(LinRange(a,b,n)))

function LogRange(a,b,n; z = 1, T = Float64)
    X = T.(range(0,1,n))
    Z = exp.(T(z) .* X)
    x = first(Z)
    y = last(Z)
    @. (T(b)-T(a))/(y-x) * (Z - x) + T(a) 
end

# previous log Range function which is false but may be still interesting
# function LogRange(a,b,n; z = 1, T = Float64)
#     X = T.(range(a,b,n))
#     a = first(X)
#     b = last(X)
#     Z = zero(X)
#     for i ∈ firstindex(X):lastindex(X)-1
#         x = X[i]
#         tmp = (b-a)/(b-x)
#         tmp = T(z) * T(log(tmp))
#         Z[i] = b - (b-a)/(tmp+1)
#     end
#     Z[end] = b
#     Z
# end

function logmesh(a, b, n; z = 1, T = Float64)
    mesh(LogRange(a, b, n; z = z, T = T))
end

@inline Base.eltype(::OneDMesh{T}) where T = T
@inline Base.eachindex(m::OneDMesh) = eachindex(m.points)
@inline Base.firstindex(m::OneDMesh) = firstindex(m.points)
@inline Base.lastindex(m::OneDMesh) = lastindex(m.points)
@inline Base.getindex(m::OneDMesh, n::Int) = m.points[n]
@inline Base.getindex(m::OneDMesh, ur::UnitRange{Int64}) = m.points[ur]
@inline Base.setindex!(m::OneDMesh{T}, val::T, n::Int) where T = m.points[n] = val
@inline Base.first(m::OneDMesh) = m[firstindex(m)]
@inline Base.last(m::OneDMesh) = m[lastindex(m)]

@inline Base.length(m::OneDMesh) = length(m.points)
@inline Base.size(m::OneDMesh) = size(m.points)

@inline points(m::OneDMesh) = m.points

@inline function findindex(m::OneDMesh, x)
    for i in eachindex(m)
        if m[i] > x
            return i-1
        end
    end
    if x == m[end]
        return lastindex(m)
    else
        return lastindex(m)+1
    end
end
   
@inline Base.iterate(m::OneDMesh, state = 1) = state > length(m) ? nothing : (m[state],state+1)
@inline left(m::OneDMesh) = m[1]
@inline right(m::OneDMesh) = m[end]
@inline edges(m::OneDMesh) = (m[1],m[end])