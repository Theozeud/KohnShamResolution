abstract type AbstractPolynomial{T} end
@inline Base.eltype(::AbstractPolynomial{T}) where T = T