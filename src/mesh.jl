abstract type AbstractMesh end

struct OneDMesh{T} <: AbstractMesh
    points::Vector{T}
    steps::Vector{T}
end

function mesh(points::AbstractArray{T})  where T<: Real
    _points = union(sort(points))
    steps = [_points[i+1] - _points[i] for i in eachindex(_points)[1:end-1]]
    OneDMesh(_points , steps)
end

linmesh(a, b, n) = mesh(LinRange(a,b,n))

function logmesh(a,b,n)
    X = LinRange(a,b,n)
    a = first(X)
    b = last(X)
    Z = zero(X)
    for i ∈ firstindex(X):lastindex(X)-1
        x = X[i]
        tmp = (b-a)/(b-x)
        tmp = log(tmp)
        Z[i] = b - (b-a)/(tmp+1)
    end
    Z[end] = b
    mesh(Z)
end

@inline Base.eltype(::OneDMesh{T}) where T = T
@inline Base.eachindex(m::OneDMesh) = eachindex(m.points)
@inline Base.firstindex(m::OneDMesh) = firstindex(m.points)
@inline Base.lastindex(m::OneDMesh) = lastindex(m.points)
@inline Base.getindex(m::OneDMesh, n::Int) = m.points[n]
@inline Base.setindex!(m::OneDMesh{T}, val::T, n::Int) where T = m.points[n] = val
@inline Base.first(m::OneDMesh) = m[firstindex(m)]
@inline Base.last(m::OneDMesh) = m[lastindex(m)]

@inline Base.length(m::OneDMesh) = length(m.points)
@inline Base.size(m::OneDMesh) = size(m.points)

@inline points(m::OneDMesh) = m.points
@inline steps(m::OneDMesh) = m.steps


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