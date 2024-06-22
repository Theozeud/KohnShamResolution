using KohnShamResolution
using Plots

# Test on integrated Legendre basis
T = Float64
Rmin = 0
cutting_pre = 50
Rmax = (1.5  + cutting_pre*log(10))
Nmesh = 5
m = linmesh(Rmin, Rmax, Nmesh)

order = 4
basis = IntLegendreBasis(m, T; order = order, left = false, right = false)

X = LinRange(Rmin, Rmax,100*Nmesh)
plt = plot(legend = :outertopright)
for i ∈ eachindex(basis)
    plot!(X,basis[i].(X))
end
plt

deriv_basis = deriv(basis)
A   = mass_matrix(deriv_basis, Rmin, Rmax)
M₀  = mass_matrix(basis, Rmin, Rmax)
M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

z = 1
l = 0
H = 1/2 * (A + l*(l+1)*M₋₂) - z .* M₋₁

using LinearAlgebra
values, vectors = eigen(H, M₀)


function check_symmetry(M)
    for I ∈ CartesianIndices(M)
        if M[I[1],I[2]] ≠ M[I[2],I[1]]
            print(I)
        end
    end
end

check_symmetry(M₀)

function mass_matrix2(lpb::LaurentPolynomialBasis, a::Real, b::Real)
    T = eltype(first(lpb))
    A = zeros(T, (length(lpb), length(lpb)))
    for I ∈ lpb.cross_index
        @inbounds A[I[1],I[2]] = scalar_product(lpb[I[1]], lpb[I[2]], a, b)
        @inbounds A[I[2],I[1]] = A[I[1],I[2]] 
    end
    A
end

M₀2 = mass_matrix(basis, Rmin, Rmax)



## Reorder legendre basis
using Memoize
@memoize function IntLegendre2(mesh::OneDMesh, i::Int, n::Int, T::Type = Float64)
    @assert i ≤ lastindex(mesh) - 1
    pₙ = Legendre(mesh[i], mesh[i+1], n, T, true)
    int_pₙ = integrate(pₙ)
    int_pₙ - int_pₙ(mesh[i])
end

function IntLegendre_element2(mesh::OneDMesh, i::Int, n::Int, T::Type = Float64)
    pₙ = IntLegendre2(mesh, i, n, T)
    PiecewiseLaurentPolynomial(mesh, [pₙ], [i], T(0))
end


function IntLegendreBasis2(mesh::OneDMesh{TM}, T::Type = Float64; order::Int = 1, left::Bool = true, right::Bool = left) where TM
    @assert order ≥ 1
    basis = PiecewiseLaurentPolynomial{T,TM}[]
    if left       
        push!(basis, HatFunctionP1(mesh, firstindex(mesh), T))
    end
    for i ∈ eachindex(mesh)[begin+1:end-1]
        push!(basis, HatFunctionP1(mesh, i, T))
    end
    if right       
        push!(basis, HatFunctionP1(mesh, lastindex(mesh), T))
    end
    for i ∈ eachindex(mesh)[begin:end-1]
        for p ∈ 2:order 
            push!(basis, IntLegendre_element2(mesh, i, p-1, T))
        end
    end
    LaurentPolynomialBasis(basis)
end

order = 3
basis2 = IntLegendreBasis2(m, T; order = order, left = false, right = false)
M2  = mass_matrix(basis2, Rmin, Rmax)