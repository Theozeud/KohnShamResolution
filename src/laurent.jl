"""
LaurentPolynomial

"""
mutable struct LaurentPolynomial
    coeffs
    degmin::Int
    haslog::Bool
    coeff_log
end

@inline deg(p::LaurentPolynomial) = (p.degmin, p.degmin+length(p.coeffs)-1)
@inline degpos(p::LaurentPolynomial) = p.degmin+length(p.coeffs)-1
@inline degneg(p::LaurentPolynomial) = p.degmin
@inline haslog(p::LaurentPolynomial) = p.haslog
@inline convert_index(p::LaurentPolynomial, i::Int) = i-degneg(p)+1
@inline ismonomial(p::LaurentPolynomial) = degneg(p) == degpos(p) && !haslog(p)
@inline Base.eachindex(p::LaurentPolynomial) = degneg(p):degpos(p)
@inline Base.getindex(p::LaurentPolynomial, i::Int) =  i ∈ eachindex(p) ? p.coeffs[i-degneg(p)+1] : 0  
@inline Base.firstindex(p::LaurentPolynomial) = degneg(p)
@inline Base.lastindex(p::LaurentPolynomial) = degpos(p)

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
    p.coeffs = p.coeffs[dn-degneg(p)+1:dp-degneg(p)+1]
    p.degmin = dn
    if p.haslog && p.coeff_log == 0
        p.haslog = false
    end
    p
end

function Base.show(io::IO, p::LaurentPolynomial)
    elag!(p)
    str = ""
    if degneg(p) ≤ -1
        for i in degneg(p):-1
            if p[i] ≠ 0
                if p[i] ≠ 1
                    str *= string(p[i])*string(" X^")*string(i)*" + "
                else
                    str *= string(" X^(")*string(i)*") + "
                end
            end
        end
    end
    if degpos(p)≥0
        if p[0] != 0 || ismonomial(p)
            str *= string(p[0])*" + "
        end
    end
    if degpos(p)≥1
        if p[1] ≠ 0
            if p[1] ≠ 1
                str *= string(p[1])*string(" X")*" + "
            else
                str *= string(" X")*" + "
            end
        end
    end
    if degpos(p)≥2
        for i in 2:degpos(p)
            if p[i] ≠ 0
                if p[i] ≠ 1
                    str *= string(p[i])*string(" X^")*string(i)*" + "
                else
                    str *= string(" X^")*string(i)*" + "
                end
            end
        end
    end
    str = str[1:length(str)-3]
    if str == "" 
        str *= "0"
    end
    println(io, str)
end


function (p::LaurentPolynomial)(x)
    y = p[end]
    for i ∈ degpos(p)-1:-1:0
        y = y*x + p[i]
    end
    z = p[begin]/x
    for i ∈ degneg(p)+1:-1
        z = (z + p[i])/x
    end
    z+y+p.coeff_log * log(x)
end



function shift!(p::LaurentPolynomial, n::Int)
    p.degmin += n
end

function Base.:*(r::Real, p::LaurentPolynomial)
    elag!(LaurentPolynomial(p.coeffs .* r, p.degmin,haslog(p), p.coeff_log))
end


function Base.:+(p::LaurentPolynomial, q::LaurentPolynomial)
    if degneg(p) ≤ degneg(q)
        if degpos(p) < degneg(q)
            new_coeffs = zeros(degpos(q) - degneg(p) + 1)
            for (i,e) in zip(1:length(p.coeffs),p)
                new_coeffs[i] = e
            end
            for (i,e) in zip(1:length(q.coeffs),q)
                new_coeffs[end - i + 1] = e
            end
            LaurentPolynomial(new_coeffs, degmin(q), haslog(p) || haslog(q), p.coeff_log + q.coeff_log)
        else
            new_coeffs = zeros(degpos(q) - degneg(p) + 1)
            for (i,e) in zip(1:length(p.coeffs),p)
                new_coeffs[i] = e
            end
            for (i,e) in zip(1:length(q.coeffs),q)
                new_coeffs[end - i + 1] = e
            end
            LaurentPolynomial(new_coeffs, degmin(q), haslog(p) || haslog(q), p.coeff_log + q.coeff_log)
        end
    else
        q+p
    end
end

#=
function Base.:*(p::Polynomial, q::Polynomial)
    elag!(p)
    elag!(q)
    new_coeffs = zeros(deg(p)+deg(q)+1)
    for i ∈ eachindex(new_coeffs)
        for j ∈ 0:i-1
            new_coeffs[i] += p[j]*q[i-j-1]
        end
    end
    Polynomial(new_coeffs)
end
=#



function integrate!(p::LaurentPolynomial)
    if haslog(p)
        @error "We can't integrate a laurent polynomial with already a log term."
    end
    if p[-1] != 0
        p.haslog = true
        p.coeff_log = p[-1]
    end
    p.coeffs = p.coeffs .* [i == -1 ? 0 : 1/1+i for i in eachindex(p)]
    shift!(p,1)
end


function integrate(p::LaurentPolynomial)
    if haslog(p)
        @error "We can't integrate a laurent polynomial with already a log term."
    end
    _haslog = p[-1] != 0 ? true : false
    coeff_log = p[-1]
    new_coeffs = p.coeffs .* [i == -1 ? 0 : 1/1+i for i in eachindex(p)]
    LaurentPolynomial(new_coeffs, p.degmin+1, _haslog, coeff_log)
end

function integrate(p::LaurentPolynomial, a::Real, b::Real)
    int_p = integrate(p)
    int_p(b) - int_p(a)
end

#=
function deriv!(p::Polynomial)
    p.coeffs = p.coeffs .* [i for i in 1:deg(p)]
    pop_deg!(p,1)    
end

function deriv(p::Polynomial)
    new_coeffs = p.coeffs[2:end] .* [i for i in 1:deg(p)]
    Polynomial(new_coeffs)
end
=#


# Test

p = LaurentPolynomial([4,5,0,1,2],-3,false,0)
q = integrate(p)
r = integrate(q)