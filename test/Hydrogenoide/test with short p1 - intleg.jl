using KohnShamResolution
using LinearAlgebra
using Plots

# Creation of the model
z = 1
N = 1
KM = KohnShamExtended(z = z,N = N)

# Choice of the method
method = ConstantODA(1.0)

# Discretization 
Nmesh = 30
lₕ = 0
Rmin = 0.00000001
cutting_pre = 10
Rmax = 300 #(1.5 * log(z) + cutting_pre*log(10))/z
m = logmesh(Rmin, Rmax, Nmesh; z= 0.5)

T = Float64
normalize = true

# Block P1
left = false
right = false
p1 = ShortP1Basis(m, T; normalize = normalize, left = left, right = right)

# Block IntLegendre
ordermin = 2
ordermax = 3
intleg = ShortIntLegendreBasis(m, T; normalize = normalize, ordermin = ordermin, ordermax = ordermax)

# Combinaison
basis = CombineShortPolynomialBasis(p1, intleg)

# Solve
deriv_basis = deriv(basis)
A = mass_matrix(deriv_basis)
M₀ = mass_matrix(basis)

#M₋₂ = weight_mass_matrix(basis, -2)


#=
# Plots of basis
X = LinRange(Rmin,Rmax,1000)
plt = plot()
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(X, p.(X))
end
plt
=#

function weight_mass_matrix_2(cb::CombineShortPolynomialBasis, weight::LaurentPolynomial)
    T = KohnShamResolution.bottom_type(first(cb))
    A = zeros(T, (length(cb), length(cb)))
    for b ∈ cb.blocks
        @views ABlock = A[KohnShamResolution.getrangerow(b), KohnShamResolution.getrangecolumn(b)]
        if KohnShamResolution.isdiagonal(b)
            KohnShamResolution.fill_weight_mass_matrix!(KohnShamResolution.getbasis(cb, KohnShamResolution.getindex(b,1)), weight, ABlock)
        else
            fill_weight_mass_matrix_2!(KohnShamResolution.getbasis(cb, KohnShamResolution.getindex(b,1)), KohnShamResolution.getbasis(cb, KohnShamResolution.getindex(b,2)), weight, b.interaction_index, A)
            @views ABlockT = A[KohnShamResolution.getrangecolumn(b), KohnShamResolution.getrangerow(b)]
            @. ABlockT = ABlock'
        end
    end
    A 
end

function weight_mass_matrix_2(cb::CombineShortPolynomialBasis, n::Int)
    weight_mass_matrix_2(cb, Monomial(n))
end


function fill_weight_mass_matrix_2!(spb1::ShortPolynomialBasis, spb2::ShortPolynomialBasis, weight::LaurentPolynomial, interaction_index::Vector{CartesianIndex{2}}, A)
    for I ∈ interaction_index
        for (i,j) ∈ KohnShamResolution.intersection_with_indices(KohnShamResolution.getsegments(spb1, I[1]), KohnShamResolution.getsegments(spb2, I[2]))
            P = KohnShamResolution.getpolynomial(spb1, I[1], i)
            Q = KohnShamResolution.getpolynomial(spb2, I[2], j)
            invϕ = KohnShamResolution.getinvshift(spb1, I[1], i)
            dinvϕ = invϕ[1]
            weight_shift = weight ∘ invϕ
            try
                @inbounds A[I[1], I[2]] += dinvϕ * KohnShamResolution.weight_scalar_product(P, Q, weight_shift, spb1.elements.binf, spb1.elements.bsup)
            catch
                @show I
                @show KohnShamResolution.getsegments(spb1, I[1])
                @show KohnShamResolution.getsegments(spb1, I[2])
                @show (i,j)
                @show P
                @show Q
                @show invϕ
                @show weight
                @show weight_shift
            end
            
        end
        if KohnShamResolution.isnormalized(spb1)
            @inbounds A[I[1], I[2]] *= KohnShamResolution.getnormalization(spb1, I[1]) 
        end
        if isnormalized(spb2)
            @inbounds A[I[1], I[2]] *= KohnShamResolution.getnormalization(spb2, I[2])
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
end


M₋₁ = weight_mass_matrix_2(basis, -1)


H = 1/2 * (A ) - z .* M₋₁ #+ lₕ*(lₕ+1)*M₋₂
ϵ, U = eigen(H,M₀)
scatter(ϵ)