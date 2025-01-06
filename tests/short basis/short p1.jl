using KohnShamResolution
using Test
using Plots

# Test with P1
T = Float64
Rmin = 0
Rmax = 5
Nmesh = 6
m = linmesh(Rmin, Rmax, Nmesh)
normalize = true
left = false
right = false

p1 = P1Basis(m, T; left = left, right = right)
@time "Complete P1 Mass matrix" mm = mass_matrix(p1)
display(mm)

shortp1 = ShortP1Basis(m, T; normalize = normalize, left = left, right = right)
@time "Short P1 Mass matrix" short_mm = mass_matrix(shortp1)
display(short_mm)

@time "Complete P1 weight Mass -1 matrix" M₋₁ = weight_mass_matrix(p1, -1)
display(M₋₁)
@time "Short P1 weight Mass -1 matrix" shortM₋₁ = weight_mass_matrix(shortp1, -1)
display(shortM₋₁)

@time "Complete P1 weight Mass -2 matrix" M₋₂ = weight_mass_matrix(p1, -2)
display(M₋₂)
@time "Short P1 weight Mass -2 matrix" shortM₋₂ = weight_mass_matrix(shortp1, -2)
display(shortM₋₂)

@time "Complete P1 deriv Mass matrix" dm = mass_matrix(deriv(p1))
display(dm)
@time "Short P1 deriv Mass matrix" shortdm = mass_matrix(deriv(shortp1))
display(shortdm)

# Plot Basis
X = LinRange(Rmin, Rmax, Nmesh * 100)
plt_basis = plot(legend =false)
for i ∈ 1:length(shortp1)
    p = build_basis(shortp1, i)
    plot!(plt_basis, X, p.(X), lw = 4)
end
xlabel!("Maillage")
plt_basis

# Plot Derivative of the Basis
deriv_basis = deriv(shortp1)
X = LinRange(Rmin, Rmax, Nmesh * 100)
plt_derivbasis = plot(legend =false)
for i ∈ 1:length(deriv_basis)
    p = build_basis(deriv_basis, i)
    plot!(plt_derivbasis, X, p.(X))
end
plt_derivbasis