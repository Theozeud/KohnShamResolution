########################################################################################
#                                  Integrated Legendre Basis
########################################################################################

struct IntLegendreBasis{TL <: LaurentPolynomial} <: AbstractLaurentPolynomialBasis
    mesh::OneDMesh
    polynomial::Vector{TL}
    order::Int
    left::Bool
    right::Bool

    function IntLegendreBasis(mesh::OneDMesh, T::Type = Float64; order::Int = 1, left::Bool = true, right::Bool = left)
        @assert order ≥ 1
        basis = LaurentPolynomial{T}[]
    
        hf1 = LaurentPolynomial([sqrt(T(3))/sqrt(T(2))], 1, false, T(0))
        hf2 = LaurentPolynomial([-sqrt(T(3))/sqrt(T(2))], 1, false, T(0))
        push!(basis, hf1)
        push!(basis, hf2)
    
        for n ∈ 2:order
            Qₙ = intLegendre(n-1; T = T, normalize = true)
            push!(basis, Qₙ)
        end
    
        new{eltype(basis)}(mesh, basis, order, left, right)
    end
    
end

@inline nb_hat_functions(ilb::IntLegendreBasis) = length(ilb.mesh) - 2 + ilb.left + ilb.right
@inline getpolynomial(ilb::IntLegendreBasis, n::Int) = ilb.polynomial[n]
@inline bottom_type(ilb::IntLegendreBasis) = eltype(first(ilb.polynomial))
@inline size(ilb::IntLegendreBasis) = (ilb.order-1) * (length(ilb.mesh) - 1) + nb_hat_functions(ilb)
@inline Order(ilb::IntLegendreBasis) = ilb.order


"""
    mass_matrix(ilb::IntLegendreBasis)

Compute the mass matrix for the Integrated-Legendre Polynomials basis. 
    
The upper-right block is the mass matrix of the P1 basis. The lower-left block is a diagonal
made of the interaction of polynomials of order bigger than one, in ascending order of 
orders, arranged mesh by mesh. The anti-diagonal block are the interaction between the P1 
basis and the other elements of the basis.
"""
function mass_matrix(ilb::IntLegendreBasis)
    @unpack left, right, mesh, polynomial, order = ilb
    nbhf = nb_hat_functions(ilb)
    siz  = size(ilb)
    T = bottom_type(ilb)
    A = zeros(T, (size(ilb), size(ilb)))

    # Block of hat function
    if left
        A[1,1] = (mesh[2] -  mesh[1])/3 * 3/2
    else
        A[1,1] = (mesh[3] -  mesh[1])/3 * 3/2
    end
    for I ∈ 2:nbhf-1
        A[I,I]   = (mesh[I+1 + !left] -  mesh[I-1 + !left])/3 * 3/2
        A[I,I-1] = (mesh[I+1 + !left] -  mesh[I + !left])/6   * 3/2
        A[I-1,I] = A[I,I-1]
    end
    A[nbhf, nbhf-1] = (mesh[nbhf + !left] -  mesh[nbhf - 1 + !left])/6 * 3/2
    A[nbhf-1, nbhf] = A[nbhf, nbhf-1]
    if right
        A[nbhf, nbhf] = (mesh[nbhf + !left] -  mesh[nbhf - 1 + !left])/3 * 3/2
    else
        A[nbhf, nbhf] = (mesh[nbhf + 1 + !left] -  mesh[nbhf + !left])/3 * 3/2
    end

    # Diagonal for polynomial of higer order
    idxmsh= 1
    icount = 1
    for I ∈ nbhf + 1 : siz
        A[I,I] = 2/(mesh[idxmsh + 1]-mesh[idxmsh])
        icount += 1
        if icount == order
            idxmsh += 1
            icount = 1
        end
    end

    # Interaction hat functions and polynomial of higher order
    hf_up = getpolynomial(ilb, 1)
    hf_down = getpolynomial(ilb, 2)

    for i ∈ eachindex(mesh)[1:end-1]
        if i == 1
            if left
                for n ∈ 2:order
                    Qₙ = getpolynomial(ilb, n+1)
                    I = nbhf + n - 1
                    A[I, i]   = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_down, Qₙ, -1, 1)
                    A[I, i+1] = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_up, Qₙ, -1 ,1)
                    A[i, I]   = A[I, i]
                    A[i+1, I] = A[I, i+1]
                end
            else
                for n ∈ 2:order
                    Qₙ = getpolynomial(ilb, n+1)
                    I = nbhf + n - 1
                    A[I, i]   = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_up, Qₙ, -1, 1)
                    A[i, I]   = A[I, i]
                end
            end
        elseif i == lastindex(mesh)-1
            for n ∈ 2:order
                Qₙ = getpolynomial(ilb, n+1)
                I = (siz - order + 1) + n - 1
                A[I, i] = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_down, Qₙ, -1 ,1)
                A[i, I]   = A[I, i]
            end
        else 
            for n ∈ 2:order
                Qₙ = getpolynomial(ilb, n+1)
                I = nbhf  + (order - 1) * (i-1) + n - 1
                A[I, i]   = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_down, Qₙ, -1, 1)
                A[I, i+1] = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_up, Qₙ, -1 ,1)
                A[i, I]   = A[I, i]
                A[i+1, I] = A[I, i+1]
            end
        end
    end
    A
end

function weight_mass_matrix(::IntLegendreBasis, ::LaurentPolynomial)
    @error "weight_mass_matrix is not implemented for integrated legendre polynomial basis."
end

function weight_mass_matrix(ilb::IntLegendreBasis, p::Int)
    @unpack left, right, mesh, polynomial, order = ilb
    nbhf = nb_hat_functions(ilb)
    T = bottom_type(ilb)
    A = zeros(T, (size(ilb), size(ilb)))

    # Block of P1 Basis
    @views AP1 = A[1:nbhf, 1:nbhf]
    hf_up = getpolynomial(ilb, 1)
    hf_down = getpolynomial(ilb, 2)
    for I ∈ 1:nbhf
        if I == 1 && !left
            shit_weight_up = shift_power_weight(T, -1, 1, mesh[1], mesh[2], p)
            shit_weight_down = shift_power_weight(T, -1, 1, mesh[2], mesh[3], p)
            AP1[I,I] = weight_scalar_product(hf_up, hf_up, shit_weight_up, -1, 1) + weight_scalar_product(hf_down, hf_down, shit_weight_down, -1, 1)
            AP1[I,I+1] = weight_scalar_product(hf_down, hf_up, shit_weight_down, -1, 1)
        elseif I == nbhf && !right
            shit_weight_up = shift_power_weight(T, -1, 1, mesh[end-2], mesh[end-1], p)
            shit_weight_down = shift_power_weight(T, -1, 1, mesh[end-1], mesh[end], p)
            AP1[I,I] = weight_scalar_product(hf_up, hf_up, shit_weight_up, -1, 1) + weight_scalar_product(hf_down, hf_down, shit_weight_down, -1, 1)
            AP1[I,I+1] = AP1[I,I-1]
        else 
            shit_weight_up = shift_power_weight(T, -1, 1, mesh[I - 1 + !left], mesh[I + !left], p)
            shit_weight_down = shift_power_weight(T, -1, 1, mesh[I + !left], mesh[I + 1 + !left], p)
            AP1[I,I] = weight_scalar_product(hf_up, hf_up, shit_weight_up, -1, 1) + weight_scalar_product(hf_down, hf_down, shit_weight_down, -1, 1)
            AP1[I,I+1] = weight_scalar_product(hf_down, hf_up, shit_weight_up, -1, 1)
            AP1[I,I-1] = AP1[I-1,I]
        end
    end

    # Block of integrated legendre polynomials
    @views AHO = A[nbhf+1:end, nbhf+1:end]
    for i ∈ eachindex(mesh)[1:end-1]
        shit_weight = shift_power_weight(T, -1, 1, mesh[i], mesh[i+1], p)
        idxm = i*(order-1)
        for n ∈ 2:order
            Qₙ = getpolynomial(ilb, n+1)
            IN = idxm + n - 1
            for m ∈ 2:order
                Qₘ = getpolynomial(ilb, m+1)
                IM = idxm + m - 1
                if m > n
                    AHO[IN, IM] = AHO[IM, IN]
                else
                    AHO[IN, IM] = weight_scalar_product(Qₙ, Qₘ, shit_weight, -1, 1)
                end
            end
        end
    end

    # Block of interaction P1 - integrated legendre polynomials
    @views AP1HO = A[nbfh+1:end, 1:nbhf]
    for i ∈ eachindex(mesh)[1:end-1]
        shift_weight = shift_power_weight(T, -1, 1, mesh[i + 1], mesh[i], p)
        idxm    = i*(order-1)
        I = i - !left
        for n ∈ 2 :order
            IN  = idxm  + n - 1
            Qₙ = getpolynomial(ilb, n+1)
            if i == firstindex(mesh) && !left
                AP1HO[IN, I+1] = weight_scalar_product(hf_up, Qₙ, shift_weight, -1, 1)
            elseif i == lastindex(mesh)-1 && !right
                AP1HO[IN, I] = weight_scalar_product(hf_up, Qₙ, shift_weight, -1, 1)
            else
                AP1HO[IN, I]   = weight_scalar_product(hf_down, Qₙ, shift_weight, -1, 1)
                AP1HO[IN, I+1] = weight_scalar_product(hf_up,   Qₙ, shift_weight, -1, 1)
            end
        end
    end
    @views sym_AP1H0 = A[1:nbhf, nbfh+1:end]
    @. sym_AP1H0 = AP1HO'

    A
end

@memoize function fast_monom_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, n::Int, a::Real, b::Real)
    scalar_product(Monomial(n) * p, q, a, b)
end

@inline function weight_scalar_product(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}, weight::LaurentPolynomial{TW}, a::Real, b::Real) where {TP, TQ, TW}
    NewT = promote_type(TP,TQ, TW)
    sum = NewT(0)
    for i ∈ eachindex(weight)
        sum += weight[i] * fast_monom_scalar_product(p, q, i, a, b)
    end
    sum
end

@inline function shift_power_weight(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real, n::Int)
    @assert n≥0
    if n == 0
        return LaurentPolynomial([T(1)], 0, false, T(0))
    else
        c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
        c0 = -T(a) * T(c1) + T(mᵢ)
        w = LaurentPolynomial([c0, c1], 0, false, 0)
        return w^n
    end
end

# Deriv
function deriv(lpb::IntLegendreBasis)
    LegendreBasis(lpb.mesh, bottom_type(lpb); order = Order(lpb)-1, lpb.left , lpb.right)
end

#= legendre polynomial rescaled on [m[i], m[i+1]]
@memoize function IntLegendre(mesh::OneDMesh, i::Int, n::Int, T::Type = Float64)
    @assert i ≤ lastindex(mesh) - 1
    pₙ = Legendre(mesh[i], mesh[i+1], n, T, false)
    int_pₙ = integrate(pₙ)
    int_pₙ - int_pₙ(mesh[i])
end

function IntLegendre_element(mesh::OneDMesh, i::Int, n::Int, T::Type = Float64)
    pₙ = IntLegendre(mesh, i, n, T)
    PiecewiseLaurentPolynomial(mesh, [pₙ], [i], T(0))
end

function IntLegendreBasis(mesh::OneDMesh{TM}, T::Type = Float64; order::Int = 1, left::Bool = true, right::Bool = left) where TM
    @assert order ≥ 1
    basis = PiecewiseLaurentPolynomial{T,TM}[]
    for i ∈ eachindex(mesh)[begin:end-1]
        if i != firstindex(mesh) || left       
            push!(basis, HatFunctionP1(mesh, i, T))
        end
        for p ∈ 2:order 
            push!(basis, IntLegendre_element(mesh, i, p, T))
        end
    end
    if right
        push!(basis, HatFunctionP1(mesh, lastindex(mesh), T))
    end
    LaurentPolynomialBasis(basis)
end
=#

#=
# Build LaurentPolynomial on basis
function build_on_basis(basis::LaurentPolynomialBasis{TL}, coeff) where TL
    if isempty(basis)
        @error "The basis is empty."
    end
    @assert length(coeff) == length(basis)
    poly = zero(first(basis.elements))
    for i ∈ eachindex(basis)
        poly += coeff[i] * basis[i]
    end
    poly
end
=#


########################################################################################
#                                  Legendre Basis
########################################################################################

struct LegendreBasis{TL <: LaurentPolynomial} <: AbstractLaurentPolynomialBasis
    mesh::OneDMesh
    polynomial::Vector{TL}
    order::Int
    left::Bool
    right::Bool

    function LegendreBasis(mesh::OneDMesh, T::Type = Float64; order::Int = 1, left::Bool = true, right::Bool = left)
        @assert order ≥ 1
        basis = LaurentPolynomial{T}[]
    
        dhf1 = LaurentPolynomial([sqrt(T(3))/sqrt(T(2))], 0, false, T(0))
        dhf2 = LaurentPolynomial([-sqrt(T(3))/sqrt(T(2))], 0, false, T(0))
        push!(basis, dhf1)
        push!(basis, dhf2)
    
        for n ∈ 2:order
            Pₙ = Legendre(n-1; T = T, normalize = true)
            push!(basis, Pₙ)
        end
    
        new{eltype(basis)}(mesh, basis, order, left, right)
    end
    
end


function mass_matrix(ilb::LegendreBasis)
    @unpack left, right, mesh, polynomial, order = ilb
    nbhf = nb_hat_functions(ilb)
    T = bottom_type(ilb)
    A = zeros(T, (size(ilb), size(ilb)))

    # Block of P1 Basis
    @views AP1 = A[1:nbhf, 1:nbhf]
    hf_up = getpolynomial(ilb, 1)
    hf_down = getpolynomial(ilb, 2)
    for I ∈ 1:nbhf
        if I == 1 && !left
            shit_weight_up = shift_power_weight(T, -1, 1, mesh[1], mesh[2], p)
            shit_weight_down = shift_power_weight(T, -1, 1, mesh[2], mesh[3], p)
            AP1[I,I] = weight_scalar_product(hf_up, hf_up, shit_weight_up, -1, 1) + weight_scalar_product(hf_down, hf_down, shit_weight_down, -1, 1)
            AP1[I,I+1] = weight_scalar_product(hf_down, hf_up, shit_weight_down, -1, 1)
        elseif I == nbhf && !right
            shit_weight_up = shift_power_weight(T, -1, 1, mesh[end-2], mesh[end-1], p)
            shit_weight_down = shift_power_weight(T, -1, 1, mesh[end-1], mesh[end], p)
            AP1[I,I] = weight_scalar_product(hf_up, hf_up, shit_weight_up, -1, 1) + weight_scalar_product(hf_down, hf_down, shit_weight_down, -1, 1)
            AP1[I,I+1] = AP1[I,I-1]
        else 
            shit_weight_up = shift_power_weight(T, -1, 1, mesh[I - 1 + !left], mesh[I + !left], p)
            shit_weight_down = shift_power_weight(T, -1, 1, mesh[I + !left], mesh[I + 1 + !left], p)
            AP1[I,I] = weight_scalar_product(hf_up, hf_up, shit_weight_up, -1, 1) + weight_scalar_product(hf_down, hf_down, shit_weight_down, -1, 1)
            AP1[I,I+1] = weight_scalar_product(hf_down, hf_up, shit_weight_up, -1, 1)
            AP1[I,I-1] = AP1[I-1,I]
        end
    end

    # Block of integrated legendre polynomials
    @views AHO = A[nbhf+1:end, nbhf+1:end]
    for i ∈ eachindex(mesh)[1:end-1]
        shit_weight = shift_power_weight(T, -1, 1, mesh[i], mesh[i+1], p)
        idxm = i*(order-1)
        for n ∈ 2:order
            Qₙ = getpolynomial(ilb, n+1)
            IN = idxm + n - 1
            for m ∈ 2:order
                Qₘ = getpolynomial(ilb, m+1)
                IM = idxm + m - 1
                if m > n
                    AHO[IN, IM] = AHO[IM, IN]
                else
                    AHO[IN, IM] = weight_scalar_product(Qₙ, Qₘ, shit_weight, -1, 1)
                end
            end
        end
    end

end
