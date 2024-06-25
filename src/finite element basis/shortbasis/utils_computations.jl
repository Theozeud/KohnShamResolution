@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    @assert n≥0
    if n == 0
        return LaurentPolynomial([T(1)], 0, false, T(0))
    else
        c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
        c0 = -T(a) * T(c1) + T(mᵢ)
        return LaurentPolynomial([c0, c1], 0, false, 0)
    end
end

@memoize function fast_monom_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, n::Int, a::Real, b::Real)
    scalar_product(Monomial(n) * p, q, a, b)
end

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TP, TQ, TW}
    NewT = promote_type(TP,TQ, TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        if i ≥ 0
            sum += weight[i] * fast_monom_scalar_product(p, q, i, a, b)
        elseif i == -1

        elseif i == -2

        else
            @error "We cant"
        end
    end
    sum
end
