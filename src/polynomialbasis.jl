struct PolynomialBasis 
    polynomials
end

@inline Base.size(pb::PolynomialBasis) = size(pb.polynomials)
@inline Base.getindex(pb::PolynomialBasis, n::Int) =  pb.polynomials[n] 
@inline Base.eachindex(pb::PolynomialBasis) = eachindex(pb.polynomials)


scalar_product(p::Polynomial, q::Polynomial) = integrate(p*q)
scalar_product(p::Polynomial, q::Polynomial, a::Real, b::Real) = integrate(p*q,a,b)


# Mass matrix
function mass_matrix(pb::PolynomialBasis, a::Real, b::Real)
    A = zeros(size(pb),size(pb))
    for i in eachindex(pb)
        for j in eachindex(pb)
            A[i,j] = scalar_product(pb[i], pb[j], a, b)
        end
    end
    A
end

