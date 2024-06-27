using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
m = linmesh(1,5,5)
normalize = true
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

@time "Complete P1 weight Mass -1 matrix" M₋₁ = weight_mass_matrix(p1, -1)
display(M₋₁)
@time "Short P1 weight Mass -1 matrix" shortM₋₁ = weight_mass_matrix(shortp1, -1)
display(shortM₋₁)

@time "Complete P1 weight Mass -2 matrix" M₋₂ = weight_mass_matrix(p1, -2)
display(M₋₂)
@time "Short P1 weight Mass -2 matrix" shortM₋₂ = weight_mass_matrix(shortp1, -2)
display(shortM₋₂)

