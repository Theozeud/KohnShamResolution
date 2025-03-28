mutable struct LaurentPolynomial{T} <:AbstractPolynomial{T}
    coeffs::Vector{T}
    degmin::Int
    haslog::Bool
    coeff_log::T
end

Monomial(deg::Int, coeff = 1) = LaurentPolynomial([coeff], deg, false, oftype(coeff, 0))
Polynomial(coeff, degmin::Int = 0) = LaurentPolynomial(coeff, degmin, false, eltype(coeff)(0))
RootsPolynomial(roots::AbstractVector) = length(roots) == 1 ?  Monomial(1) - first(roots) : (Monomial(1) - first(roots)) * RootsPolynomial(roots[2:end])

RandMonomial(deg::Int) = Monomial(deg)
RandMonomial(T::Type, deg::Int) = Monomial(deg, rand(T))
RandPolynomial(degmax::Int, degmin::Int = 0) = Polynomial(rand(degmax - degmin +1), degmin)
RandPolynomial(T::Type, degmax::Int, degmin::Int = 0) = Polynomial(rand(T, degmax - degmin +1), degmin)

@inline _convert(::Type{T}, p::LaurentPolynomial) where T = LaurentPolynomial(T.(p.coeffs), p.degmin, p.haslog, T(p.coeff_log))

@inline Base.eachindex(p::LaurentPolynomial) = degmin(p):degmax(p)
@inline Base.getindex(p::LaurentPolynomial, i::Int) =  i ∈ eachindex(p) ? p.coeffs[i-degmin(p)+1] : 0
@inline Base.getindex(p::LaurentPolynomial, ur::UnitRange{Int64}) = p.points[ur .+ 1 .- degmin(p)]
@inline Base.setindex!(p::LaurentPolynomial{T}, val::T, i::Int) where T =  p.coeffs[i-degmin(p)+1] = val   
@inline Base.firstindex(p::LaurentPolynomial) = degmin(p)
@inline Base.lastindex(p::LaurentPolynomial) = degmax(p)

@inline deg(p::LaurentPolynomial) = (p.degmin, p.degmin + length(p.coeffs) - 1)
@inline degmax(p::LaurentPolynomial) = p.degmin + length(p.coeffs) - 1
@inline degmin(p::LaurentPolynomial) = p.degmin
@inline haslog(p::LaurentPolynomial) = p.haslog
@inline ismonomial(p::LaurentPolynomial) = degmin(p) == degmax(p) && !haslog(p)

@inline Base.zero(p::LaurentPolynomial{T}) where T = LaurentPolynomial(zero(p.coeffs), p.degmin, false, T(0))
@inline Base.zero(::Type{LaurentPolynomial{T}}, args...) where T = Laurent_zero(T, 0, 0)
function Laurent_zero(T::Type, degmin::Int, degmax::Int)
    @assert degmax ≥ degmin
    LaurentPolynomial(zeros(T, degmax-degmin+1), degmin, false, T(0))
end
function Base.iszero(p::LaurentPolynomial)
    elag!(p)
    p.coeffs == [0] && p.degmin == 0 && !haslog(p)
end

Base.round(p::LaurentPolynomial; digits = 14) = LaurentPolynomial(round.(p.coeffs; digits = digits), degmin(p), haslog(p), p.coeff_log)

function elag!(p::LaurentPolynomial)
    (dn,dp) = deg(p)
    while p[dp] == 0
        if dp == dn
            break
        end
        dp -= 1
    end
    while p[dn] == 0
        if dn == dp
            break
        end
        dn += 1
    end
    x = p[dn]
    p.coeffs = p.coeffs[dn-degmin(p)+1:dp-degmin(p)+1]
    if p.haslog && iszero(p.coeff_log) 
        p.haslog = false
    end
    if dn == dp && iszero(x)
        p.degmin = 0
    else
        p.degmin = dn
    end
    p
end

function Base.show(io::IO, p::LaurentPolynomial)
    elag!(p)
    str = ""
    if degmin(p) ≤ -1
        for i in degmin(p):-1
            if p[i] ≠ 0
                if p[i] ≠ 1
                    str *= string(p[i])*" X^"*string(i)*" + "
                else
                    str *= "X^("*string(i)*") + "
                end
            end
        end
    end
    if degmax(p)≥0
        if p[0] != 0
            str *= string(p[0])*" + "
        end
    end
    if degmax(p)≥1
        if p[1] ≠ 0
            if p[1] ≠ 1
                str *= string(p[1])*" X + "
            else
                str *= "X + "
            end
        end
    end
    if degmax(p)≥2
        for i in 2:degmax(p)
            if p[i] ≠ 0
                if p[i] ≠ 1
                    str *= string(p[i])*" X^"*string(i)*" + "
                else
                    str *= "X^"*string(i)*" + "
                end
            end
        end
    end

    if haslog(p)
        str *= string(p.coeff_log)*" log(X)"
    else
        str = str[1:length(str)-3]
    end

    if str == "" 
        str *= "0"
    end
    print(io, str)
end

function (p::LaurentPolynomial)(x)
    y = 0
    if degmax(p) ≥ 0
        y = p[end]
        for i ∈ degmax(p)-1:-1:0
            y = y*x + p[i]
        end
    end
    z = 0
    invx = inv(x)
    if degmin(p) < 0
        z = p[begin] * invx
        for i ∈ degmin(p)+1:-1
            z = (z + p[i]) * invx
        end
    end
    if haslog(p)
        z+y+p.coeff_log * log(x)
    else
        z+y
    end
end

function shift!(p::LaurentPolynomial, n::Int)
    p.degmin += n
end

function shift(p::LaurentPolynomial, n::Int)
    LaurentPolynomial(p.coeffs, p.degmin + n, haslog(p), p.coeff_log)
end

##################################################################################
#                            Elementary Computations
##################################################################################

function Base.:-( p::LaurentPolynomial)
    coeffs = - p.coeffs
    coeff_log = - p.coeff_log
    LaurentPolynomial(coeffs, p.degmin, haslog(p), coeff_log)
end

function Base.:*(r::T, p::LaurentPolynomial{TP}) where {T, TP}
    NewT = promote_type(T, TP)
    if iszero(r)
        return LaurentPolynomial(NewT[0], 0, false, NewT(0))
    elseif r == 1
        return _convert(NewT,p)
    else
        coeffs = NewT.(p.coeffs) .* NewT(r)
        coeff_log = NewT(p.coeff_log) * NewT(r)
        return LaurentPolynomial(coeffs, p.degmin, haslog(p), coeff_log)
    end
end

function Base.:*(p::LaurentPolynomial, r::Real)
    r * p
end

function Base.:/(p::LaurentPolynomial, r::Real)
    if iszero(r)
        throw(DivideError())
    end
    return p * (1 /r)
end


function Base.:+(p::LaurentPolynomial{TP}, x::T) where {TP, T}
    NewT = promote_type(TP, T)
    degmin_r = min(degmin(p), 0)
    degmax_r = max(degmax(p), 0)
    r = Laurent_zero(NewT, degmin_r, degmax_r)
    for i in eachindex(r)
        if i == 0
            r[i] = NewT(p[i]) + NewT(x)
        else
            r[i] = NewT(p[i])
        end
    end
    r.coeff_log = NewT(p.coeff_log)
    r.haslog = haslog(p)
    elag!(r)
    return r
end

function Base.:+(x, p::LaurentPolynomial) 
    p + x
end

function Base.:-(x, p::LaurentPolynomial) 
    x + (-p)
end

function Base.:-(p::LaurentPolynomial, x)
    p + (-x)
end

function Base.:-(p::LaurentPolynomial, q::LaurentPolynomial)
    p + (-q)
end


function Base.:+(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where {TP, TQ}
    NewT = promote_type(TP, TQ)
    degmin_r = min(degmin(p), degmin(q))
    degmax_r = max(degmax(p), degmax(q))
    r = Laurent_zero(NewT, degmin_r, degmax_r)
    @inbounds for i ∈ eachindex(r)
        r[i] = NewT(p[i]) + NewT(q[i])
    end
    r.coeff_log = NewT(p.coeff_log + q.coeff_log)
    r.haslog = p.haslog || q.haslog
    elag!(r)
    return r
end


function Base.:*(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where{TP, TQ}
    if haslog(p) || haslog(q)
        @error "We can't multiply two laurent polynomial if at least one of them have a log term."
    end
    NewT = promote_type(TP,TQ)
    r = Laurent_zero(NewT, degmin(p) + degmin(q), degmax(p)+ degmax(q))
    for i ∈ eachindex(r)
        for j ∈ eachindex(p)
            r[i] += p[j]*q[i-j]
        end
    end
    elag!(r)
end

function Base.:^(p::LaurentPolynomial{T}, n::Int) where T
    if haslog(p)
        throw(ArgumentError("Cannot power Laurent polynomials with logarithmic terms."))
    end
    if n == 0
        return LaurentPolynomial([T(1)], 0, false, T(0))
    elseif n == 1
        return p
    else
        q = p^(n÷2)
        if iseven(n) 
            return q * q
        else
            return p * q * q
        end
    end
end

function diveucl(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where{TP, TQ}
    @assert !haslog(p) && !haslog(q)
    @assert degmin(p) ≥ 0 && degmin(q) ≥ 0
    @assert !iszero(q)
    NewT = promote_type(TP,TQ)
    if degmax(p) < degmax(q)
        return (Monomial(0, NewT(0)), p)
    end
    _q = q/q[end]
    _q = _q /_q[end] #To enforce the dominant coefficient to be one if it is 0.99999999...
    Q = Monomial(0, NewT(0))
    R = _convert(NewT, p)
    while degmax(R) ≥ degmax(q)
        Q += Monomial(degmax(R) - degmax(q), NewT(R[end]))
        R -= _q * Monomial(degmax(R) - degmax(q), R[end])
    end
    return (round(Q/NewT(q[end])), R)
end


##################################################################################
#                            Integration & Derivation
##################################################################################

function integrate(p::LaurentPolynomial{T}) where T
    if haslog(p)
        @error "We can't integrate a laurent polynomial with already a log term."
    end
    _haslog = p[-1] != 0 ? true : false
    coeff_log = p[-1]
    new_coeffs = p.coeffs .* [i == -1 ? 0//1 : 1//(1+i) for i in eachindex(p)]
    LaurentPolynomial(new_coeffs, p.degmin+1, _haslog, eltype(new_coeffs)(coeff_log))
end


function integrate(p::LaurentPolynomial, a::Real, b::Real)
    int_p = integrate(p)
    int_p(b) - int_p(a)
end

function deriv(p::LaurentPolynomial{T}) where T
    new_coeffs = p.coeffs .* [i for i in eachindex(p)]
    if haslog(p)
        new_coeffs[-degmin(p)+1] = p.coeff_log
    end
    LaurentPolynomial(new_coeffs, p.degmin-1, false, T(0))
end


scalar_product(p::LaurentPolynomial, q::LaurentPolynomial) = integrate(p*q)
scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, a::Real, b::Real) = integrate(p*q, a, b)

scalar_product(p::LaurentPolynomial, x::Real, a::Real, b::Real) = integrate(p, a, b) * x
scalar_product(x::Real, p::LaurentPolynomial, a::Real, b::Real) = x * integrate(p, a, b)
scalar_product(y::Real, x::Real, a::Real, b::Real) = x*y*(b-a)

function normL2(p::LaurentPolynomial{T}, a::Real, b::Real) where T
    if iszero(p)
        return T(0)
    else
        int = scalar_product(p, p, a, b)
        if int > 0
            return sqrt(int)
        else
            return sqrt(-int)
        end
    end
end