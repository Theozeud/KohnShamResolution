abstract type AbstractMesh end

struct 1DMesh{TP,TS} <: AbstractMesh
    points
    steps
end

function mesh(points::AbstractArray{T} where T<: Real)
    _points = union(sort(points))
    steps = [_points[i+1] -  _points[i] for i in eachindex(_points)[1:end-1]]
    1Dmesh(_points , steps)
end

function mesh(point::Real, fun::Base.Callable, N::Int)
    points = zeros(N)
    points[1] = point
    for i in 2:N
        @inbounds points[i] = fun(points[i-1])
    end
    mesh(points)
end

@inline Base.length(m::1DMesh) = length(m.points)
@inline Base.size(m::1DMesh) = size(m.points)
@inline points(m::1DMesh) = m.points
@inline steps(m::1DMesh) = m.steps
@inline Base.getindex(m::1DMesh, n::Int) = m.points[n]
@inline Base.iterate(m::1DMesh, state = 1) = state > size(m) ? nothing : (m[state],state+1)
@inline left(m::1DMesh) = m[1]
@inline right(m::1DMesh) = m[end]
@inline edges(m::1DMesh) = (m[1],m[end])