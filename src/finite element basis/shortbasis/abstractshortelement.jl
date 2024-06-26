########################################################################################
#                                  Abstract Short Elements
########################################################################################
abstract type AbstractShortElements{N, T} end

@inline Base.eltype(::AbstractShortElements{N, T}) where T = T
@inline isnormalized(::AbstractShortElements{N, T}) where {N,T} = N
@inline Base.length(elem::AbstractShortElements) = elem.size
@inline getpolynomial(elem::AbstractShortElements, n::Int) = elem[n]
@inline Base.first(elem::AbstractShortElements) = elem[firstindex(elem)]