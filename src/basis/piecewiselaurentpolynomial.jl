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
@inline getfunction(pwlp::PiecewiseLaurentPolynomial, i::Int) = i ∈ pwlp.index ? pwlp.functions[findfirst(item->item == i, pwlp.index)] : x->pwlp.default_value

(pwlp::PiecewiseLaurentPolynomial)(x) = pwlp[KohnShamResolution.findindex(pwlp.mesh, x)][2](x) 


## Functions

function Base.:+(p::PiecewiseLaurentPolynomial{TP}, q::PiecewiseLaurentPolynomial{TQ}) where {TP,TQ}
    if p.mesh.points != q.mesh.points
        @error "We can't add for the moment piecewise laurent polynomial on different meshes." 
    end
    NewT = promote_type(TP,TQ)
    index_p = Set(p.index)
    index_q = Set(q.index)
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    for i ∈ index_p
        fp = getfunction(p,i)
        if i ∈ index_q
            fq = getfunction(q,i)
            sumfpfq = fp + fq
            KohnShamResolution.elag!(sumfpfq)
            if !iszero(sumfpfq)
                push!(laurent_poly, sumfpfq)
                push!(index,i)
            end
        else
            push!(laurent_poly, fp)
            push!(index,i)
        end
    end
    for i ∈ setdiff(index_q,index_p)
        fq = getfunction(q,i)
        push!(laurent_poly, fq)
        push!(index,i)
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly,index, NewT(p.default_value) + NewT(q.default_value))
end

function Base.:*(p::PiecewiseLaurentPolynomial{TP}, q::PiecewiseLaurentPolynomial{TQ}) where{TP, TQ}
    if p.mesh.points != q.mesh.points
        @error "We can't multiply for the moment piecewise laurent polynomial on different meshes." 
    end
    NewT = promote_type(TP,TQ)
    index_p = Set(p.index)
    index_q = Set(q.index)
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    for i ∈ index_p
        fp = getfunction(p,i)
        if i ∈ index_q
            fq = getfunction(q,i)
            prodfpfq = fp * fq
            KohnShamResolution.elag!(prodfpfq)
            if !iszero(sumfpfq)
                push!(laurent_poly, prodfpfq)
                push!(index,i)
            end
        else
            push!(laurent_poly, fp)
            push!(index,i)
        end
    end
    for i ∈ setdiff(index_q,index_p)
        fq = getfunction(q,i)
        push!(laurent_poly * p.default_value, fq)
        push!(index,i)
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly,index, NewT(p.default_value) * NewT(q.default_value))
end

function integrate!(p::PiecewiseLaurentPolynomial)
    for i ∈ eachindex(p)
        integrate!(p.functions[findfirst(item->item == i, p.index)])
    end
end


function integrate(p::PiecewiseLaurentPolynomial{T}) where T
    if haslog(p)
        @error "We can't integrate a laurent polynomial with already a log term."
    end
    _haslog = p[-1] != 0 ? true : false
    coeff_log = p[-1]
    new_coeffs = p.coeffs .* [i == -1 ? 0//1 : 1//(1+i) for i in eachindex(p)]
    LaurentPolynomial(new_coeffs, p.degmin+1, _haslog, eltype(new_coeffs)(coeff_log))
end

function integrate(p::PiecewiseLaurentPolynomial, a::Real, b::Real)
    int_p = integrate(p)
    int_p(b) - int_p(a)
end

function deriv!(p::PiecewiseLaurentPolynomial)
    p.coeffs = p.coeffs .* [i for i in eachindex(p)]
    shift!(p,-1)
    if haslog(p)
        p.coeffs[-degmin(p)] = p.coeff_log
        p.haslog = false
    end
    p
end

function deriv(p::PiecewiseLaurentPolynomial)
    new_coeffs = p.coeffs .* [i for i in eachindex(p)]
    if haslog(p)
        new_coeffs[-degmin(p)+1] = p.coeff_log
    end
    LaurentPolynomial(new_coeffs, p.degmin-1, false, 0.0)
end


scalar_product(p::PiecewiseLaurentPolynomial, q::PiecewiseLaurentPolynomial) = integrate(p*q)
scalar_product(p::PiecewiseLaurentPolynomial, q::PiecewiseLaurentPolynomial, a::Real, b::Real) = integrate(p*q,a,b)
