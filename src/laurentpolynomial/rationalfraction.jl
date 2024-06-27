mutable struct RationalFraction{T} <: AbstractPolynomial{T}
    ent::LaurentPolynomial{T}
    num::LaurentPolynomial{T}
    denom::LaurentPolynomial{T}
end

function Base.:/(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where {TP, TQ}
    @assert !haslog(p) && !haslog(q)
    NewT = promote_type(TP,TQ)
    if ismonomial(q)
        return shift(p, -degmin(q))/q[begin]
    end
    degm = min(degmin(p), degmin(q))
    _p = shift(p , -degm)
    _q = shift(q , -degm)
    Q,R = diveucl(_p, _q)
    RationalFraction(Q, R/NewT(q[end]), convert(NewT, q)/NewT(q[end]))
end

@inline Base.eltype(::RationalFraction{T}) where T = T

@inline deg_ent(rf::RationalFraction) = deg(rf.ent)
@inline degmax_ent(rf::RationalFraction) = degmax(rf.ent)
@inline degmin_ent(rf::RationalFraction) = degmin(rf.ent)
@inline deg_num(rf::RationalFraction) = deg(rf.num)
@inline degmax_num(rf::RationalFraction) = degmax(rf.num)
@inline degmin_num(rf::RationalFraction) = degmin(rf.num)
@inline deg_denom(rf::RationalFraction) = deg(rf.denom)
@inline degmax_denom(rf::RationalFraction) = degmax(rf.denom)
@inline degmin_denom(rf::RationalFraction) = degmin(rf.denom)

function (rf::RationalFraction)(x)
    rf.ent(x) + rf.num(x)/rf.denom(x)
end


##################################################################################
#                            Elementary Computations
##################################################################################

#=
function elag!(rf::RationalFraction{T})
    elag!(rf.num)
    if iszero(rf.num)
        @warn "Numerator is null but the object is still a rational fraction."
    end
    elag!(rf.denom)
    Q, R = diveucl(rf.num, rf.denom)
    if degmax(R) == 0
        rf.ent += Q + rf.num/R[0]
        elag!(ent)
        rf.num = zero(rf.num)
    else
        rf.denom = R
        rf.ent += Q
        elag!(ent)
    end
end
=#

function Base.:+(rf::RationalFraction{TR}, p::LaurentPolynomial{TP}) where {TR, TP}
    NewT = promote_type(TR,TP)
    rf.ent + p
    RationalFraction(rf.ent + p, convert(NewT, rf.num), convert(NewT, rf.denom))
end

function Base.:+(p::LaurentPolynomial, rf::RationalFraction)
    rf + p
end

function Base.:*(rf::RationalFraction, p::LaurentPolynomial)
    (p * rf.num) / rf.denom + rf.ent * p
end

function Base.:*(p::LaurentPolynomial, rf::RationalFraction)
    rf * p
end

##################################################################################
#                            Integration & Derivation
##################################################################################

function integrate(rf::RationalFraction, a::Real, b::Real)
    if degmax_denom(rf) ≥ 3
        @error "No analatycal expression for integrating rational fractions with a denominator of degree higher than 2."
    elseif degmax_denom(rf) == 1
        B = rf.num[0]
        D = rf.denom[1]
        E = rf.denom[0]
        @assert -E/D > b || -E/D < a 
        return integrate(rf.ent, a, b) + B/D * log((D*b + E)/(D*a + E))
    elseif degmax_denom(rf) == 2
        A = rf.num[1]
        B = rf.num[0]
        C = rf.denom[2]
        D = rf.denom[1]
        E = rf.denom[0]
        Δ = D^2 - 4*E*C
        if Δ > 0
            r₁ = (-D - sqrt(Δ))/(2*C)
            r₂ = (-D + sqrt(Δ))/(2*C)
            @assert (r₁ > b) || (r₂ < a) || (r₂ > b && r₁ < a)
        elseif Δ == 0
            r₀ = -D/(2*C)
            @assert (r₀ > b) || (r₀ < a)
        end
        C1 = D/(2*C)
        C2 = (4*E*C - D^2)/(4*C^2)
        C3 = B - A*C1
        if C2 > 0
            sqrtC2 = sqrt(C2)
            return A/(2*C) * log( abs( ((b+C1)^2 + C2) / ((a+C1)^2 + C2) ) ) + C3/(C*sqrtC2) * (atan((b + C1)/sqrtC2) - atan((a + C1)/sqrtC2))
        else
            sqrtC2 = sqrt(-C2)
            return A/(2*C) * log( abs( ((b+C1)^2 + C2) / ((a+C1)^2 + C2) ) ) + C3/(C*2*sqrtC2) * (log(abs((B-a)/(A-a))) - log(abs((B+a)/(A+a))))
        end
    end
end

