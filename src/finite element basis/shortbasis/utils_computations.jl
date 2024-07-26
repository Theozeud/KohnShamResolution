@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
    c0 = -T(a) * T(c1) + T(mᵢ)
    LaurentPolynomial([c0, c1], 0, false, T(0))
end

############################################################################
# Fast Scalar Product
############################################################################

@memoize function fast_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, a::Real, b::Real)
    scalar_product(p, q, a, b)
end

@memoize function fast_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, a::Real, b::Real)
    scalar_product(p, q*l, a, b)
end

@memoize function fast_monom_scalar_product(p::LaurentPolynomial, n::Int, a::Real, b::Real)
    scalar_product(Monomial(n), p, a, b)
end

@memoize function fast_monom_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, n::Int, a::Real, b::Real)
    scalar_product(Monomial(n) * p, q, a, b)
end

@memoize function fast_monom_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, n::Int, a::Real, b::Real)
    scalar_product(Monomial(n)*p, q*l, a, b)
end

############################################################################
# Weight Scalar Product
############################################################################

# Polynomial

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, weight::LaurentPolynomial{TW},a::Real, b::Real) where {TP, TW}
    NewT = promote_type(TP,TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * fast_monom_scalar_product(p, i, a, b)
    end
    sum
end
@inline weight_scalar_product(p::LaurentPolynomial, weight::LaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(p, weight∘ϕ, a, b)

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TP, TQ, TW}
    NewT = promote_type(TP,TQ, TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * fast_monom_scalar_product(p, q, i, a, b)
    end
    sum
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::LaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, weight∘ϕ, a, b)

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, l::LaurentPolynomial{TL}, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TP, TQ, TL, TW}
    NewT = promote_type(TP,TQ,TL,TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * fast_monom_scalar_product(p, q, l, i, a, b)
    end
    sum
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, weight::LaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, l, weight∘ϕ, a, b)


# Piecewisepolynomial

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, weight::PiecewiseLaurentPolynomial{TW}, a::Real, b::Real) where {TP, TW}
    NewT = promote_type(TP,TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        if i ∈ weight.index
            sum += weight_scalar_product(p, weight[i], a, b)
        elseif !iszero(weight.default_value)
            sum += fast_monom_scalar_product(p, 0, a, b) * weight.default_value
        end
    end
    sum
end
@inline weight_scalar_product(p::LaurentPolynomial, weight::PiecewiseLaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(p, x->(weight∘ϕ)(x), a, b)

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::PiecewiseLaurentPolynomial{TW}, a::Real, b::Real) where {TP, TQ, TW}
    NewT = promote_type(TP,TQ,TW)
    sum = NewT(0)
    sum += integrate(p*q*weight, a, b)
    #=
    for i ∈ eachindex(weight)
        if i ∈ weight.index
            sum += weight_scalar_product(p, q, weight[i], a, b)
        elseif !iszero(weight.default_value)
            sum += fast_monom_scalar_product(p, q, 0, a, b) * weight.default_value
        end
    end
    =#
    sum
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::PiecewiseLaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, x->(weight∘ϕ)(x), a, b)

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, l::LaurentPolynomial{TL}, weight::PiecewiseLaurentPolynomial{TW}, a::Real, b::Real) where {TP, TQ, TL, TW}
    NewT = promote_type(TP,TQ,TL,TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        if i ∈ weight.index
            sum += weight_scalar_product(p, q, l, weight[i], a, b)
        elseif !iszero(weight.default_value)
            sum += fast_monom_scalar_product(p, q, l, 0, a, b) * weight.default_value
        end
    end
    sum
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, weight::PiecewiseLaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, l, x->(weight∘ϕ)(x), a, b)

# RationalFraction

@inline function weight_scalar_product(p::LaurentPolynomial, weight::RationalFraction, a::Real, b::Real)
    if degmax(weight.denom) > 2
        return integrate((p*weight.num)/weight.denom, a, b; geomfun = false, enforceNullDelta = false)
    else
        return integrate((p*weight.num)/weight.denom, a, b; geomfun = true, enforceNullDelta = true)
    end
end
@inline weight_scalar_product(p::LaurentPolynomial, weight::RationalFraction, a::Real, b::Real, ϕ) = weight_scalar_product(p, weight∘ϕ, a, b)

@inline function weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::RationalFraction, a::Real, b::Real)
    if degmax(weight.denom) > 2
        return integrate((p*q*weight.num)/weight.denom, a, b; geomfun = false, enforceNullDelta = false)
    else
        return integrate((p*q*weight.num)/weight.denom, a, b; geomfun = true, enforceNullDelta = true)
    end
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::RationalFraction, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, weight∘ϕ, a, b)

@inline function weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, weight::RationalFraction, a::Real, b::Real)
    if degmax(weight.denom) > 2
        return integrate((p*q*l*weight.num)/weight.denom, a, b; geomfun = false, enforceNullDelta = false)
    else
        return integrate((p*q*l*weight.num)/weight.denom, a, b; geomfun = true, enforceNullDelta = true)
    end
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, weight::RationalFraction, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, l, weight∘ϕ, a, b)

# SommeRationalFraction

@inline function weight_scalar_product(p::LaurentPolynomial, weight::SommeRationalFraction, a::Real, b::Real)
    val = weight_scalar_product(p, weight.ent, a, b)
    for i ∈ eachindex(weight)
        val += weight_scalar_product(p, weight[i], a, b)
    end
    val
end
@inline weight_scalar_product(p::LaurentPolynomial, weight::SommeRationalFraction, a::Real, b::Real, ϕ) = weight_scalar_product(p, weight∘ϕ, a, b)

@inline function weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::SommeRationalFraction, a::Real, b::Real)
    val = weight_scalar_product(p, q, weight.ent, a, b)
    for i ∈ eachindex(weight)
        val += weight_scalar_product(p, q, weight[i], a, b)
    end
    val
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::SommeRationalFraction, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, weight∘ϕ, a, b)

@inline function weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, weight::SommeRationalFraction, a::Real, b::Real)
    val = weight_scalar_product(p, q, l, weight.ent, a, b)
    for i ∈ eachindex(weight)
        val += weight_scalar_product(p, q, l, weight[i], a, b)
    end
    val
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, weight::SommeRationalFraction, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, l, weight∘ϕ, a, b)

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, weight::Base.Callable, a::Real, b::Real) where TP
    f(x) = weight(x) * p(x)
    approximate_integral(f, (a,b); method = QuadGKJL(), reltol =  100*eps(TP), abstol =  100*eps(TP))
end
@inline weight_scalar_product(p::LaurentPolynomial, weight::Base.Callable, a::Real, b::Real, ϕ) = weight_scalar_product(p, weight∘ϕ, a, b)

# Callable

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::Base.Callable, a::Real, b::Real) where {TP, TQ}
    f(x) = weight(x) * p(x) * q(x)
    NewT = promote_type(TP,TQ)
    approximate_integral(f, (a,b); method = QuadGKJL(), reltol = 100*eps(NewT), abstol =  100*eps(NewT))
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, weight::Base.Callable, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, weight∘ϕ, a, b)

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, l::LaurentPolynomial{TL}, weight::Base.Callable, a::Real, b::Real) where {TP, TQ, TL}
    f(x) = weight(x) * p(x) * q(x) * l(x)
    NewT = promote_type(TP,TQ,TL)
    approximate_integral(f, (a,b); method = QuadGKJL(), reltol = 100*eps(NewT), abstol =  100*eps(NewT))
end
@inline weight_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, weight::Base.Callable, a::Real, b::Real, ϕ) = weight_scalar_product(p, q, l, weight∘ϕ, a, b)