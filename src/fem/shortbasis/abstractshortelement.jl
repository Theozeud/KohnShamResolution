########################################################################################
#                                  Abstract Short Elements
########################################################################################
abstract type AbstractShortElements{T} end

@inline Base.eltype(::AbstractShortElements{T}) where T = T
@inline Base.length(elem::AbstractShortElements) = elem.size
@inline getpolynomial(elem::AbstractShortElements, n::Int) = elem[n]
@inline getderivpolynomial(elem::AbstractShortElements, n::Int) = getderivpolynomial(elem)[n]
@inline Base.first(elem::AbstractShortElements) = elem[firstindex(elem)]