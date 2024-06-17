# Using the Bonnet Recurrence formula to compute the nth, normalized or not, Legendre Polynomial scaled on [a,b]
@memoize function Legendre(a::Real, b::Real, n::Int, T::Type = Float64, normalize = false)
    if n == 0
        if normalize
            return Polynomial([T(1)], 0) * sqrt(T(1)/T(b-a))
        else
            Polynomial([T(1)], 0)
        end
    elseif n == 1
        if normalize
            return Polynomial([T(-b-a)/T(b-a), T(2)/T(b-a)]) * sqrt(T(3)/T(b-a))
        else
            return Polynomial([T(-b-a)/T(b-a), T(2)/T(b-a)])
        end
    else
        coeff = normalize ? sqrt(T(2*n + 1)/T(b-a)) : T(1)
        pₙ₋₁ = Legendre(a, b, n-1, T, false)
        pₙ₋₂ = Legendre(a, b, n-2, T, false)
        pₙ₋₁ = (2*n-1)/n * Polynomial([T(-b-a)/T(b-a), T(2)/T(b-a)]) * pₙ₋₁ - (n-1)/n * pₙ₋₂
        return pₙ₋₁ * coeff
    end
end
