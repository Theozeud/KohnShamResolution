abstract type AbstractMesh end

struct OneDMesh{T} <: AbstractMesh
    points::Vector{T}
end

function mesh(points::AbstractArray{T})  where T<: Real
    _points = union(sort(points))
    OneDMesh(_points)
end

linmesh(a, b, n; T = Float64) = mesh(T.(LinRange(a,b,n)))

function LogRange(a,b,n; z = 1, T = Float64)
    X = T.(range(0,1,n))
    Z = exp.(T(z) .* X)
    x = first(Z)
    y = last(Z)
    @. (T(b)-T(a))/(y-x) * (Z - x) + T(a) 
end

function logmesh(a, b, n; z = 1, T = Float64)
    mesh(LogRange(a, b, n; z = z, T = T))
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

# Mesh adapt to x->Cx*e^(-x)
# For that, one have to solve y = Cxe^(-x), i.e
# (-x)e^(-x) = -y/C with y < C/e
# Using th Lambert W function, this gives
# x = -W(-y/C, -1) ou x = -W(-y/C, 0)
function LinearExpRange(a, b, n; T = Float64, coeff = one(T))
    @assert a < b
    ta = Base.convert(T, a)
    tb = Base.convert(T, b)
    eval_a = coeff * ta * exp(-ta)
    eval_b = coeff * tb * exp(-tb)
    if ta ≤ tb ≤ 1
        Y1 = T.(range(eval_a, eval_b, n))
        return lambertw.(-Y1 ./coeff, 0)
    elseif ta < 1 < tb
        eval_one = coeff*exp(-one(T))
        n1 = Int(round(n * (eval_one-eval_a)/(2*eval_one-eval_a - eval_b)))
        n2 = n - n1
        Y1 = T.(range(eval_a, eval_one, n1))
        Y2 = T.(range(eval_one, eval_b, n2 + 1))
        return reduce(vcat, (-lambertw.(-Y1 ./coeff, 0), -lambertw.(-Y2[2:end] ./coeff, -1)))
    elseif 1 ≤ a ≤ b
        Y2 = T.(range(eval_b, eval_a, n))
        return lambertw.(Y2, -1)
    end
end

function linearexpmesh(a, b, n; T = Float64, coeff = one(T))
    mesh(LinearExpRange(a, b, n; T = T, coeff = coeff))
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