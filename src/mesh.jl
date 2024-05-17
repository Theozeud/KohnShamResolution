abstract type AbstractMesh end

struct OneDMesh{TP,TS} <: AbstractMesh
    points::TP
    steps::TS
end

function mesh(points::AbstractArray{T} where T<: Real)
    _points = union(sort(points))
    steps = [_points[i+1] - _points[i] for i in eachindex(_points)[1:end-1]]
    OneDmesh(_points , steps)
end

function mesh(point::Real, fun::Base.Callable, N::Int)
    points = zeros(N)
    points[1] = point
    for i in 2:N
        @inbounds points[i] = fun(points[i-1])
    end
    mesh(points)
end

@inline Base.length(m::OneDMesh) = length(m.points)
@inline Base.size(m::OneDMesh) = size(m.points)
@inline points(m::OneDMesh) = m.points
@inline steps(m::OneDMesh) = m.steps
@inline Base.getindex(m::OneDMesh, n::Int) = m.points[n]
@inline Base.iterate(m::OneDMesh, state = 1) = state > size(m) ? nothing : (m[state],state+1)
@inline left(m::OneDMesh) = m[1]
@inline right(m::OneDMesh) = m[end]
@inline edges(m::OneDMesh) = (m[1],m[end])