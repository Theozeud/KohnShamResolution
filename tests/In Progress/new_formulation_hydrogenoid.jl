using KohnShamResolution
include("../../benchmarktools/hydrogenoid/setup.jl")


@inline function KohnShamResolution.weight_scalar_product(x::Real, y::Real, weight::PiecewiseLaurentPolynomial{TW}, a::Real, b::Real) where {TW}
    NewT = promote_type(typeof(x),typeof(y),TW)
    sum = NewT(0)
    sum += x * y * integrate(weight, a, b)
    sum
end
@inline KohnShamResolution.weight_scalar_product(x::Real, y::Real, weight::PiecewiseLaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(x, y, z->(weight∘ϕ)(z), a, b)


@inline function KohnShamResolution.weight_scalar_product(x::Real, q::LaurentPolynomial{TQ}, weight::PiecewiseLaurentPolynomial{TW}, a::Real, b::Real) where {TQ, TW}
    NewT = promote_type(typeof(x),TQ,TW)
    sum = NewT(0)
    sum += x *integrate(q*weight, a, b)
    sum
end
@inline KohnShamResolution.weight_scalar_product(x::Real, q::LaurentPolynomial, weight::PiecewiseLaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(x, q, z->(weight∘ϕ)(z), a, b)

@inline function KohnShamResolution.weight_scalar_product(p::LaurentPolynomial{TP}, y::Real, weight::PiecewiseLaurentPolynomial{TW}, a::Real, b::Real) where {TP, TW}
    NewT = promote_type(TP,typeof(y),TW)
    sum = NewT(0)
    sum += y * integrate(p*weight, a, b)
    sum
end
@inline KohnShamResolution.weight_scalar_product(p::LaurentPolynomial, y::Real, weight::PiecewiseLaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(p, y, x->(weight∘ϕ)(x), a, b)

using KohnShamResolution: fast_monom_scalar_product
@inline function KohnShamResolution.weight_scalar_product(x::Real, y::Real, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TW}
    NewT = promote_type(typeof(x),typeof(y), TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * KohnShamResolution.fast_monom_scalar_product(Polynomial([x]), Polynomial([y]), i, a, b)
    end
    sum
end
@inline KohnShamResolution.weight_scalar_product(x::Real, y::Real, weight::LaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(x, y, weight∘ϕ, a, b)

@inline function KohnShamResolution.weight_scalar_product(x::Real, q::LaurentPolynomial{TQ}, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TQ, TW}
    NewT = promote_type(typeof(x),TQ, TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * KohnShamResolution.fast_monom_scalar_product(Polynomial([x]), q, i, a, b)
    end
    sum
end
@inline KohnShamResolution.weight_scalar_product(x::Real, q::LaurentPolynomial, weight::LaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(x, q, weight∘ϕ, a, b)

@inline function KohnShamResolution.weight_scalar_product(p::LaurentPolynomial{TP}, y::Real, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TP, TW}
    NewT = promote_type(TP,typeof(y), TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * KohnShamResolution.fast_monom_scalar_product(p, Polynomial([y]), i, a, b)
    end
    sum
end
@inline KohnShamResolution.weight_scalar_product(p::LaurentPolynomial, y::Real, weight::LaurentPolynomial, a::Real, b::Real, ϕ) = weight_scalar_product(p, y, weight∘ϕ, a, b)


# Creation of missing matrix

function weight_stiffness_matrix(spb::ShortPolynomialBasis, weight)
    @unpack elements, mesh, size = spb
    T = bottom_type(spb)
    A = zeros(T, size, size)
    fill_weight_stiffness_matrix!(spb, weight, A)
    A
end

function weight_stiffness_matrix(spb::ShortPolynomialBasis, n::Int)
    weight_stiffness_matrix(spb, Monomial(n))
end

function fill_weight_stiffness_matrix!(spb::ShortPolynomialBasis, weight, A)
    @unpack elements, mesh = spb
    for I ∈ spb.coupling_index
        for (i,j) ∈ KohnShamResolution.intersection_with_indices(KohnShamResolution.getsegments(spb, I[1]), KohnShamResolution.getsegments(spb, I[2]))
            P = KohnShamResolution.getderivpolynomial(spb, I[1], i)
            Q = KohnShamResolution.getderivpolynomial(spb, I[2], j)
            ϕ = KohnShamResolution.getshift(spb, I[1], i)
            dϕ = ϕ[1]
            invϕ = KohnShamResolution.getinvshift(spb, I[1], i)
            @inbounds A[I[1], I[2]] += dϕ * KohnShamResolution.weight_scalar_product(P, Q, weight, spb.elements.binf, spb.elements.bsup, invϕ)
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end

function special_matrix(spb::ShortPolynomialBasis, weight)
    @unpack elements, mesh, size = spb
    T = bottom_type(spb)
    A = zeros(T, size, size)
    fill_special_matrix!(spb, weight, A)
    A
end

function special_matrix(spb::ShortPolynomialBasis, n::Int)
    special_matrix(spb, Monomial(n))
end

function fill_special_matrix!(spb::ShortPolynomialBasis, weight, A)
    @unpack elements, mesh = spb
    for I ∈ spb.coupling_index
        for (i,j) ∈ KohnShamResolution.intersection_with_indices(KohnShamResolution.getsegments(spb, I[1]), KohnShamResolution.getsegments(spb, I[2]))
            Pprime  = KohnShamResolution.getderivpolynomial(spb, I[1], i)
            Qprime  = KohnShamResolution.getderivpolynomial(spb, I[2], j)
            P       = KohnShamResolution.getpolynomial(spb, I[1], i)
            Q       = KohnShamResolution.getpolynomial(spb, I[2], j)
            invϕ = KohnShamResolution.getinvshift(spb, I[1], i)
            @inbounds A[I[1], I[2]] +=  KohnShamResolution.weight_scalar_product(Pprime, Q, weight, spb.elements.binf, spb.elements.bsup, invϕ) + 
            KohnShamResolution.weight_scalar_product(P, Qprime, weight, spb.elements.binf, spb.elements.bsup, invϕ)
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end


# Eigenvalues and Eigenvectors
function eigvals_hydro2(problem)
    @unpack T, z, l, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis = problem
    mesh = typemesh(zero(T), Rmax, Nmesh; T = T, optsmesh...)
    basis = typebasis(mesh, T; optsbasis...)
    A₂  = Symmetric(weight_stiffness_matrix(basis, 2))
    M₀  = Symmetric(mass_matrix(basis))
    M₁  = Symmetric(weight_mass_matrix(basis, 1))
    M₂  = Symmetric(weight_mass_matrix(basis, 2))
    G   = Symmetric(special_matrix(basis, 1))
    if l == 0
        H       = T(0.5) * (M₀ + G + A₂ - 2*z*M₁)
        Λ       =  eigvals(inv(M₂) * H)
        return HydrogenoidSolution(problem, Λ, [1])
    else
        H       = T(0.5) * ((l*(l+1)+1)*M₀ + G + A₂ - 2*z*M₁)
        Λ       =  eigvals(inv(M₂) * H)
        return HydrogenoidSolution(problem, Λ, [1])
    end
end

function convergenceNmesh2(vecNmesh::AbstractVector, problems; nums = [1])
    ΔΛ = Dict()
    ΔU = Dict()
    for problem ∈ problems
        @unpack T, name = problem
        println(name)
        Δλ = zeros(T, length(vecNmesh), length(nums))
        @inbounds for i ∈ eachindex(vecNmesh)
            newprob = HydrogenoidProblem(problem; Nmesh = vecNmesh[i], nλ = nums, nU = nums)
            @time "Nmesh = $(vecNmesh[i])" sol = eigvals_hydro2(newprob)
            Δλ[i,:] .= sol.Δλ
        end
        ΔΛ[name] = Δλ
    end
    HydrogenoidConvergenceNmesh(problems, vecNmesh, ΔΛ, ΔU, nums)
end

# CONVERGENCE WITH RESPECT TO NMESH

probgeomesh = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 60, 
                            Nmesh         = 40,
                            typemesh      = linmesh, 
                            typebasis     = ShortP1Basis, 
                            optsmesh      = (;),  
                            optsbasis     = (left = true, right = false,),                           
                            name          = "IntLeg2-geomesh",
                            nU            = nothing)
                        
error = convergenceNmesh2(2 .^(3:8), [probgeomesh])

convergence_plot_Nmesh(error)