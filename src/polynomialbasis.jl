mutable struct Polynomial 
    coeffs
end

@inline deg(p::Polynomial) = length(p.coeffs)-1

@inline Base.getindex(p::Polynomial, n::Int) =  try p.coeffs[n] 
                                                catch 
                                                    0 end

@inline Base.eachindex(p::Polynomial) = eachindex(p.coeffs)

@inline Base.lastindex(p::Polynomial) = lastindex(p.coeffs)

function Base.show(io::IO, p::Polynomial)
    elag!(p)
    str = ""
    if deg(p) == 0
        str *= string(p[1])
    elseif p[1] != 0
        str *= string(p[1]) * " + "
    end
    if deg(p) >= 1
        if p[2] ≠ 0
            if p[2] ≠ 1
                str *= string(p[2])*string(" X")
            else
                str *= string(" X")
            end
            if deg(p) ≠ 1
                str *=" + "
            end
        end
    end
    for i in 2:deg(p)-1
        if p[i+1] ≠ 0
            if p[i+1] ≠ 1
                str *= string(p[i+1])*string(" X^")*string(i)*" + "
            else
                str *= string(" X^")*string(i)*" + "
            end
        end
    end
    deg(p) > 1 ? str *= string(p[end])*string(" X^")*string(deg(p)) : nothing
    println(io, str)
end

function (p::Polynomial)(x)
    y = p[end]
    for i ∈ deg(p)-1:-1:1
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
    m = deg(p)+1
    while p[m] == 0
        if m == 0
            break
        end
        m -= 1
    end
    p.coeffs = p.coeffs[1:m]
end

function Base.:*(r::Real, p::Polynomial)
    Polynomial(p.coeffs .* r)
end

function Base.:+(p::Polynomial, q::Polynomial)
    if deg(p) >= deg(q)
        new_coeffs = zero(q.coeffs)
        for i ∈ eachindex(q)
            new_coeffs[i] = p[i] .+ q[i]
        end
        Polynomial(reduce(vcat, (new_coeffs, p.coeffs[deg(q)+2:end])))
    else
        q+p
    end
end

function Base.:*(p::Polynomial, q::Polynomial)
    Polynomial(p.coeffs .* v)
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