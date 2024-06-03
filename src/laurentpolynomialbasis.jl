struct LaurentPolynomialBasis{T}
    elements::Vector{T}
end

@inline Base.size(lpb::LaurentPolynomialBasis) = size(lpb.elements)
@inline Base.getindex(lpb::LaurentPolynomialBasis, n::Int) =  lpb.elements[n] 
@inline Base.eachindex(lpb::LaurentPolynomialBasis) = eachindex(lpb.elements)


# Mass matrix
function mass_matrix(lpb::LaurentPolynomialBasis{T})
    A = zeros(eltype(lpb.elements), size(lpb),size(lpb))
    for i in eachindex(lpb)
        for j in eachindex(lpb)
            A[i,j] = scalar_product(lpb[i], lpb[j])
        end
    end
    A
end


function mass_matrix(lpb::LaurentPolynomialBasis, a::Real, b::Real)
    A = zeros(eltype(lpb.elements), size(lpb),size(lpb))
    for i in eachindex(lpb)
        for j in eachindex(lpb)
            A[i,j] = scalar_product(lpb[i], lpb[j], a, b)
        end
    end
    A
end

