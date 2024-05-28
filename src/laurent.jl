"""
    Laurent

"""
mutable struct Laurent 
    coeffs_plus
    coeffs_minus
    haslog::Bool
    coeff_log
end

@inline deg(p::Polynomial) = length(p.coeffs)-1
@inline Base.getindex(p::Polynomial, n::Int) =  n>deg(p) ? 0 : p.coeffs[n+1] 
@inline Base.eachindex(p::Polynomial) = 0:deg(p)
@inline Base.lastindex(p::Polynomial) = lastindex(p.coeffs)-1

function Base.show(io::IO, p::Polynomial)
    elag!(p)
    str = ""
    if deg(p) == 0
        str *= string(p[0])
    elseif p[0] != 0
        str *= string(p[0]) * " + "
    end
    if deg(p) >= 1
        if p[1] ≠ 0
            if p[1] ≠ 1
                str *= string(p[1])*string(" X")
            else
                str *= string(" X")
            end
            if deg(p) ≠ 1
                str *=" + "
            end
        end
    end
    for i in 2:deg(p)-1
        if p[i] ≠ 0
            if p[i] ≠ 1
                str *= string(p[i])*string(" X^")*string(i)*" + "
            else
                str *= string(" X^")*string(i)*" + "
            end
        end
    end
    deg(p) > 1 ? str *= string(p[deg(p)])*string(" X^")*string(deg(p)) : nothing
    println(io, str)
end

function (p::Polynomial)(x)
    y = p[end]
    for i ∈ deg(p)-1:-1:0
        y = y*x + p[i]
    end
    y
end

function push_deg!(p::Polynomial, n::Int)
    prepend!(p.coeffs, zeros(n)...)
end

function pop_deg!(p::Polynomial, n::Int)
    p.coeffs = p.coeffs[1+n:deg(p)+1]
end

function elag!(p::Polynomial)
    m = deg(p)
    while p[m] == 0
        if m == -1
            break
        end
        m -= 1
    end
    p.coeffs = p.coeffs[1:m+1]
end

function Base.:*(r::Real, p::Polynomial)
    Polynomial(p.coeffs .* r)
end

function Base.:+(p::Polynomial, q::Polynomial)
    if deg(p) >= deg(q)
        new_coeffs = zero(q.coeffs)
        for i ∈ eachindex(q)
            new_coeffs[i+1] = p[i] .+ q[i]
        end
        Polynomial(reduce(vcat, (new_coeffs, p.coeffs[deg(q)+2:end])))
    else
        q+p
    end
end

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