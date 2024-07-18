struct RationalFraction{T} <: AbstractPolynomial{T}
    num::LaurentPolynomial{T}
    denom::LaurentPolynomial{T}
    function RationalFraction(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where {TP, TQ}
        @assert !iszero(q)
        NewT = promote_type(TP,TQ)
        if haslog(p) || haslog(q)
            return new{NewT}(p, q)
        else
            degm = min(degmin(p), degmin(q))
            num = convert(NewT, shift(p , -degm))
            denom = convert(NewT, shift(q , -degm))
            return new{NewT}(num, denom)
        end
    end
end

@inline Base.eltype(::RationalFraction{T}) where T = T
@inline convert(::Type{T}, rf::RationalFraction) where T  = RationalFraction(convert(T, rf.num), convert(T, rf.denom))

function (rf::RationalFraction)(x)
    rf.num(x)/rf.denom(x)
end

function Base.:/(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where {TP, TQ}
    if ismonomial(q) && !haslog(p)
        return shift(p, -degmin(q))/q[begin]
    else
        return RationalFraction(p, q)
    end
end

struct SommeRationalFraction{T} <: AbstractPolynomial{T}
    ent::LaurentPolynomial
    ratiofrac::Vector{RationalFraction}
    function SommeRationalFraction(ent, ratiofrac)
        T = promote_type(eltype(ent), eltype.(ratiofrac)...)
        new{T}(ent, ratiofrac)
    end
end

@inline Base.eltype(::SommeRationalFraction{T}) where T = T

@inline Base.length(srf::SommeRationalFraction) = length(srf.ratiofrac)
@inline Base.firstindex(srf::SommeRationalFraction) = firstindex(srf.ratiofrac)
@inline Base.lastindex(srf::SommeRationalFraction) = lastindex(srf.ratiofrac)
@inline Base.eachindex(srf::SommeRationalFraction) = eachindex(srf.ratiofrac)
@inline Base.getindex(srf::SommeRationalFraction, i::Int) = srf.ratiofrac[i]

function (srf::SommeRationalFraction)(x)
    val = srf.ent(x)
    for i ∈ eachindex(srf)
        val += srf[i](x)
    end
    val
end

function inEntpart(rf::RationalFraction)
    @assert !haslog(rf.num) && !haslog(rf.denom)
    Q,R = diveucl(rf.num, rf.denom)
    newrf = RationalFraction(R, rf.denom)
    SommeRationalFraction(Q, [newrf])
end

function inEntpart(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where {TP, TQ}
    inEntpart(RationalFraction(p, q))
end

##################################################################################
#                            Elementary Computations
##################################################################################

function fraction_decomp(rf::RationalFraction)
    @error "The partial fraction decomposition has not be coded."
end

#=
function Base.:+(rf::RationalFraction, p::LaurentPolynomial)
    SommeRationalFraction(p, [rf])
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
=#


##################################################################################
#                            Integration & Derivation
##################################################################################

function integrate(rf::RationalFraction{T}, a::Real, b::Real; geomfun = false, enforceNullDelta = false) where T
    @assert !haslog(rf.num) && !haslog(rf.denom)
    if iszero(rf.num) || (degmax(rf.num) == 0 && abs(rf.num[0]) < sqrt(eps(typeof(rf.num[0]))))
        NewT = promote_type(T, typeof(a), typeof(b))
        return zero(NewT)
    end
    simplification = iszero(diveucl(rf.num, rf.denom)[2])
    if geomfun && !simplification
        if degmax(rf.denom) > 2
            @error "Impossibility to integrate a rational fraction with geometric functions for denominator of degree more than three."
        elseif degmax(rf.denom) == 2
            NewT = promote_type(T, typeof(a), typeof(b))
            val = zero(NewT)
            for k ∈ eachindex(rf.num)
                val += rf.num[k] * _integration_monome_over_deg2(k, rf.denom[2], rf.denom[1], rf.denom[0], a, b; enforceNullDelta = enforceNullDelta)
            end 
            return val
        elseif degmax(rf.denom) == 1
            NewT = promote_type(T, typeof(a), typeof(b))
            val = zero(NewT)
            for k ∈ eachindex(rf.num)
                val += rf.num[k] * _integration_monome_over_deg1(k, rf.denom[1],rf.denom[0], a, b)
            end 
            return val
        elseif degmax(rf.denom) == 0
            return integrate(rf.num, a, b)/rf.denom[0]
        end
    else
        if degmax(rf.num) ≥ degmax(rf.denom) 
            return integrate(inEntpart(rf), a, b)
        else
            if degmax(rf.denom) > 2
                return integrate(fraction_decomp(rf), a, b)
            elseif degmax(rf.denom) == 2
                A = rf.num[1]
                B = rf.num[0]
                C = rf.denom[2]
                D = rf.denom[1]
                E = rf.denom[0]
                return _integrate_partial12(A, B, C, D, E, a, b) 
            elseif degmax(rf.denom) == 1
                A = rf.num[0]
                B = rf.denom[1]
                C = rf.denom[0]
                return _integrate_partial01(A, B, C, a, b)
            end
        end
    end
end

function integrate(srf::SommeRationalFraction, a::Real, b::Real; geomfun = false)
    val = integrate(srf.ent, a, b)
    for i ∈ eachindex(srf)
        val += integrate(srf[i], a, b; geomfun = geomfun)
    end
    val
end