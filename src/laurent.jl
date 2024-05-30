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
@inline degnull(p::LaurentPolynomial) = p.degmin≤0 ? -p.degmin+1 : nothing
@inline haslog(p::LaurentPolynomial) = p.haslog
@inline ismonomial(p::LaurentPolynomial) = length(p.coeffs) == 1 && !haslog(p)
@inline Base.eachindex(p::LaurentPolynomial) = degneg(p):degpos(p)
@inline Base.getindex(p::LaurentPolynomial, n::Int) =  n ∈ eachindex(p) ? p.coeffs[n-degneg(p)+1] : 0  
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
    z+y
end


#=
function push_deg!(p::Polynomial, n::Int)
    prepend!(p.coeffs, zeros(n)...)
end

function pop_deg!(p::Polynomial, n::Int)
    p.coeffs = p.coeffs[1+n:deg(p)+1]
end
=#


function Base.:*(r::Real, p::LaurentPolynomial)
    elag!(LaurentPolynomial(p.coeffs .* r, p.degmin,haslog(p), p.coeff_log))
end


function Base.:+(p::LaurentPolynomial, q::LaurentPolynomial)
    new_coeffs = zero(max(degpos(p),degpos(q))-min(degneg(p),degneg(q))+1)
    for i in Set(eachindex(q))∩Set(eachindex(p))
    if deg(p) >= deg(q)
        for i ∈ eachindex(q)
            new_coeffs[i+1] = p[i] .+ q[i]
        end
        LaurentPolynomial(reduce(vcat, (new_coeffs, p.coeffs[deg(q)+2:end])))
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

function integrate!(p::Polynomial)
    p.coeffs = p.coeffs .* [1/i for i in 1:(deg(p)+1)]
    push_deg!(p,1)
end

function integrate(p::Polynomial)
    new_coeffs = p.coeffs .* [1/i for i in 1:(deg(p)+1)]
    Polynomial(prepend!(new_coeffs, 0))
end

function integrate(p::Polynomial, a::Real, b::Real)
    int_p = integrate(p)
    int_p(b) - int_p(a)
end

function deriv!(p::Polynomial)
    p.coeffs = p.coeffs .* [i for i in 1:deg(p)]
    pop_deg!(p,1)    
end

function deriv(p::Polynomial)
    new_coeffs = p.coeffs[2:end] .* [i for i in 1:deg(p)]
    Polynomial(new_coeffs)
end
=#