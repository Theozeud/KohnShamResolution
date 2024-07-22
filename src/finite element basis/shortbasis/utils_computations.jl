@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
    c0 = -T(a) * T(c1) + T(mᵢ)
    LaurentPolynomial([c0, c1], 0, false, T(0))
end

@memoize function fast_monom_scalar_product(p::LaurentPolynomial, n::Int, a::Real, b::Real)
    scalar_product(Monomial(n), p, a, b)
end

@memoize function fast_monom_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, n::Int, a::Real, b::Real)
    scalar_product(Monomial(n) * p, q, a, b)
end

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TP, TW}
    @assert degmin(weight) ≥ 0
    NewT = promote_type(TP,TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * fast_monom_scalar_product(p, i, a, b)
    end
    sum
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


@inline function weight_scalar_product(p::LaurentPolynomial, weight::RationalFraction, a::Real, b::Real)
    if degmax(weight.denom) > 2
        return integrate((p*weight.num)/weight.denom, a, b; geomfun = false, enforceNullDelta = false)
    else
        return integrate((p*weight.num)/weight.denom, a, b; geomfun = true, enforceNullDelta = true)
    end
end

@inline function weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::RationalFraction, a::Real, b::Real)
    if degmax(weight.denom) > 2
        return integrate((p*q*weight.num)/weight.denom, a, b; geomfun = false, enforceNullDelta = false)
    else
        return integrate((p*q*weight.num)/weight.denom, a, b; geomfun = true, enforceNullDelta = true)
    end
end

@inline function weight_scalar_product(p::LaurentPolynomial, weight::SommeRationalFraction, a::Real, b::Real) 
    val = weight_scalar_product(p, weight.ent, a, b)
    for i ∈ eachindex(weight)
        val += weight_scalar_product(p, weight[i], a, b)
    end
    val
end

@inline function weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::SommeRationalFraction, a::Real, b::Real) 
    val = weight_scalar_product(p, q, weight.ent, a, b)
    for i ∈ eachindex(weight)
        val += weight_scalar_product(p, q, weight[i], a, b)
    end
    val
end

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, weight::Base.Callable, a::Real, b::Real) where TP
    f(x) = weight(x) * p(x)
    approximate_integral(f, (a,b); method = QuadGKJL(), reltol =  100*eps(TP), abstol =  100*eps(TP))
end

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::Base.Callable, a::Real, b::Real) where {TP, TQ}
    f(x) = weight(x) * p(x) * q(x)
    NewT = promote_type(TP,TQ)
    approximate_integral(f, (a,b); method = QuadGKJL(), reltol = 100*eps(NewT), abstol =  100*eps(NewT))
end