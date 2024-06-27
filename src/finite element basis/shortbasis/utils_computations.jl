@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
    c0 = -T(a) * T(c1) + T(mᵢ)
    LaurentPolynomial([c0, c1], 0, false, 0)
end

@memoize function fast_monom_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, n::Int, a::Real, b::Real)
    scalar_product(Monomial(n) * p, q, a, b)
end

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TP, TQ, TW}
    @assert degmin(weight) ≥ 0
    NewT = promote_type(TP,TQ, TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * fast_monom_scalar_product(p, q, i, a, b)
    end
    sum
end

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::RationalFraction{TW}, a::Real, b::Real) where {TP, TQ, TW}
    @assert degmax_denom(weight) ≤ 2
    integrate(p*q/weight.denom, a,b) + weight_scalar_product(p, q, weight.ent, a, b)
end
