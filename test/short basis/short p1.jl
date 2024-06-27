using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
m = logmesh(1,5,5)
normalize = false
left = false
right = false

p1 = P1Basis(m, T; left = left, right = right)
@time "Complete P1 Mass matrix" mm = mass_matrix(p1)
display(mm)

shortp1 = ShortP1Basis(m, T; normalize = normalize, left = left, right = right)
@time "Short P1 Mass matrix" short_mm = mass_matrix(shortp1)
display(short_mm)

X = LinRange(1,5,1000)
plt = plot()
for i ∈ 1:length(shortp1)
    p = build_basis(shortp1, i)
    plot!(X, p.(X))
end
plt


#=
function weight_mass_matrix2(spb::ShortPolynomialBasis, weight::LaurentPolynomial)
    T = KohnShamResolution.bottom_type(spb)
    A = zeros(T, spb.size, spb.size)
    fill_weight_mass_matrix2!(spb, weight, A)
    A
end

function weight_mass_matrix2(spb::ShortPolynomialBasis, n::Int)
    weight_mass_matrix2(spb, Monomial(n))
end

function fill_weight_mass_matrix2!(spb::ShortPolynomialBasis, weight::LaurentPolynomial, A)
    for I ∈ spb.coupling_index
        for (i,j) ∈ KohnShamResolution.intersection_with_indices(KohnShamResolution.getsegments(spb, I[1]), KohnShamResolution.getsegments(spb, I[2]))
            @show (i,j)
            @show P = getpolynomial(spb, I[1], i)
            @show Q = getpolynomial(spb, I[2], j)
            @show ϕ = KohnShamResolution.getshift(spb, I[1], i)
            @show weight_shift = weight ∘ ϕ
            @show dϕ = ϕ[1]
            @show @inbounds A[I[1], I[2]] += dϕ * KohnShamResolution.weight_scalar_product(P, Q, weight_shift, spb.elements.binf, spb.elements.bsup)
        end
        @inbounds A[I[1], I[2]] *= KohnShamResolution.getnormalization(spb, I[1]) * KohnShamResolution.getnormalization(spb, I[2])
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end
=#

@time "Complete P1 weight Mass -1 matrix" M₋₁ = weight_mass_matrix(p1, -1)
display(M₋₁)
@time "Short P1 weight Mass -1 matrix" shortM₋₁ = weight_mass_matrix(shortp1, -1)
display(shortM₋₁)

#@time "Complete P1 weight Mass -2 matrix" M₋₂ = weight_mass_matrix(p1, -2)
#display(M₋₂)
#@time "Short P1 weight Mass -2 matrix" shortM₋₂ = weight_mass_matrix(shortp1, -2)
#display(shortM₋₂)
