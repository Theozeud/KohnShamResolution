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


# Mass matrix
function mass_matrix(ilb::IntLegendreBasis)
    @unpack left, right, mesh, polynomial, order = ilb

    T = bottom_type(ilb)
    A = zeros(T, (size(ilb), size(ilb)))

    #Block of hat function
    if left
        A[1,1] = (mesh[2] -  mesh[1])/3
    else
        A[1,1] = (mesh[3] -  mesh[1])/3
    end
    for I ∈ 2:nb_hat_functions(ilb)-1
        A[I,I]   = (mesh[I+1 + left] -  mesh[I-1 + !left])/3
        A[I,I-1] = (mesh[I+1 + left] -  mesh[I + !left])/6
        A[I-1,I] = A[I,I-1]
    end

    # Diagonal for polynomial of higer order
    idxmsh= 1
    icount = 1
    for I ∈ nb_hat_functions(ilb)+1:s
        A[I,I] = 2/(mesh[idxmsh + 1]-mesh[idxmsh])
        icount += 1
        if icount == order
            idxmsh += 1
            icount = 1
        end
    end

    # Interaction hat functions and polynomial of higher order
    for i ∈ eachindex(mesh)[1:end-1]
        for n ∈ 2:order
            Qₙ = getpolynomial(ilb, n+1)
            hf_up = getpolynomial(ilb, 1)
            hf_down = getpolynomial(ilb, 2)
            A[nb_hat_functions(ilb)+n-1, i]   = 2/(mesh[i + 1]-mesh[i]) * integrate(hf_up, Qₙ, -1, 1)
            A[nb_hat_functions(ilb)+n-1, i+1] = 2/(mesh[i + 1]-mesh[i]) * integrate(hf_down, Qₙ, -1 ,1)
        end
    end
    A
end

#=
function weight_mass_matrix(lpb::IntLegendreBasis, weight::LaurentPolynomial, a::Real, b::Real)
    T = eltype(first(lpb))
    A = zeros(T, (length(lpb), length(lpb)))
    for I ∈ lpb.cross_index
        @inbounds A[I[1],I[2]] = scalar_product(weight * lpb[I[1]], lpb[I[2]], a, b)
        @inbounds A[I[2],I[1]] = A[I[1],I[2]] 
    end
    A
end

function weight_mass_matrix(lpb::IntLegendreBasis, n::Int, a::Real, b::Real)
    weight_mass_matrix(lpb, Monomial(n), a::Real, b::Real)
end
=#
#=
# Deriv
function deriv!(lpb::LaurentPolynomialBasis)
    for i in eachindex(lpb)
        deriv!(lpb[i])
    end
    lpb
end

function deriv(lpb::LaurentPolynomialBasis)
    TL = eltype(lpb.elements)
    deriv_laurent = TL[]
    for i in eachindex(lpb)
        push!(deriv_laurent, deriv(lpb[i]))
    end
    LaurentPolynomialBasis(deriv_laurent)
end

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



