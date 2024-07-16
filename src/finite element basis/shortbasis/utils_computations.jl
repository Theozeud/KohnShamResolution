@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
    c0 = -T(a) * T(c1) + T(mᵢ)
    LaurentPolynomial([c0, c1], 0, false, T(0))
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

@inline function integration_monome_over_linear(k, A, B, a, b)
    @assert !iszero(A) || !iszero(B)
    @assert k ≥ 0
    if iszero(B)
        if k == 0
            1/A * log(abs(b/a))
        else
            1/(A * k) * (b^k - a^k)
        end
    elseif iszero(A)
        return 1/(B * (k+1)) * (b^(k+1) - a^(k+1))
    else
        C = B/A
        right =  b^(k+1)*pFq((1,k+1),(k+2,),-b/C)
        left = a^(k+1)*pFq((1,k+1),(k+2,),-a/C)
        return (right-left)/(A*(C*k+C))
    end
end

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::RationalFraction{TW}, a::Real, b::Real) where {TP, TQ, TW}
    if degmax_denom(weight) == 1
        val = weight_scalar_product(p, q, weight.ent, a, b)
        num = p*q*weight.num
        Q,R = diveucl(num, weight.denom)
        if iszero(R)
            return val += integrate(Q, a,b)
        else
            for k ∈ eachindex(num)
                val += num[k] * integration_monome_over_linear(k, weight.denom[1], weight.denom[0], a, b)
            end 
        return val
    end
    else
        return integrate(p*q*weight.num/weight.denom, a,b) + weight_scalar_product(p, q, weight.ent, a, b)
    end
end
