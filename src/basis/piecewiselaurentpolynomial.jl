struct PiecewiseLaurentPolynomial{T,TM}
    mesh::OneDMesh{TM}
    functions::Vector{LaurentPolynomial{T}}
    index::Vector{Int}
    default_value::T
end

@inline Base.eachindex(pwlp::PiecewiseLaurentPolynomial) = eachindex(pwlp.mesh)
@inline Base.firstindex(pwlp::PiecewiseLaurentPolynomial) = firstindex(pwlp.mesh)
@inline Base.lastindex(pwlp::PiecewiseLaurentPolynomial) = lastindex(pwlp.mesh)
@inline Base.getindex(pwlp::PiecewiseLaurentPolynomial, i::Int) = i ∈ pwlp.index ? (pwlp.mesh[i], pwlp.functions[findfirst(item->item == i, pwlp.index)]) : (pwlp.mesh[i], x->pwlp.default_value) 
@inline getmesh(pwlp::PiecewiseLaurentPolynomial, i::Int) = pwlp.mesh[i]
@inline getfunction(pwlp::PiecewiseLaurentPolynomial, i::Int) = i ∈ pwlp.index ? pwlp.functions[i] : pwlp.default_value

(pwlp::PiecewiseLaurentPolynomial)(x) = pwlp[KohnShamResolution.findindex(pwlp.mesh, x)][2](x) 


## Functions

function Base.:+(p::PiecewiseLaurentPolynomial{TP}, q::PiecewiseLaurentPolynomial{TQ}) where {TP,TQ}
    NewT = promote_type(TP,TQ)
    r = Laurent_zero(NewT, min(degmin(p), degmin(q)), max(degmax(p), degmax(q)))
    for i in eachindex(r)
        r[i] = NewT(p[i]) + NewT(q[i])
    end
    r.coeff_log = p.coeff_log + q.coeff_log
    r.haslog = p.haslog || q.haslog
    r
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
    r
end

function integrate!(p::LaurentPolynomial)
    if haslog(p)
        @error "We can't integrate a laurent polynomial with already a log term."
    end
    if p[-1] != 0
        p.haslog = true
        p.coeff_log = p[-1]
    end
    p.coeffs = p.coeffs .* [i == -1 ? 0 : 1//(1+i) for i in eachindex(p)]
    shift!(p,1)
end


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

function deriv!(p::LaurentPolynomial)
    p.coeffs = p.coeffs .* [i for i in eachindex(p)]
    shift!(p,-1)
    if haslog(p)
        p.coeffs[-degmin(p)] = p.coeff_log
        p.haslog = false
    end
    p
end

function deriv(p::LaurentPolynomial)
    new_coeffs = p.coeffs .* [i for i in eachindex(p)]
    if haslog(p)
        new_coeffs[-degmin(p)+1] = p.coeff_log
    end
    LaurentPolynomial(new_coeffs, p.degmin-1, false, 0.0)
end


scalar_product(p::LaurentPolynomial, q::LaurentPolynomial) = integrate(p*q)
scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, a::Real, b::Real) = integrate(p*q,a,b)





## Hat function
function HatFunctionP1(mesh::OneDMesh, i::Int, T::Type = Float64)
    i+=1
    if i == firstindex(mesh)
        pc = mesh[i]
        pr = mesh[i+1]
        right = LaurentPolynomial([pr/(pr-pc),T(1)/(pc-pr)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [right], [i+1], T(0))
    elseif i == lastindex(mesh)
        pl = mesh[i-1]
        pc = mesh[i]
        left = LaurentPolynomial([pl/(pl-pc),T(1)/(pc-pl)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [left], [i], T(0))
    else
        pl = mesh[i-1]
        pc = mesh[i]
        pr = mesh[i+1]
        left = LaurentPolynomial([pl/(pl-pc),T(1)/(pc-pl)], 0, false, T(0))
        right = LaurentPolynomial([pr/(pr-pc),T(1)/(pc-pr)], 0, false, T(0))
        PiecewiseLaurentPolynomial(mesh, [left, right], [i,i+1], T(0))
    end
end


    