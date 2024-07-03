using KohnShamResolution
using Plots

# Parameters of the discretization
T = Float64
Rmin = 0
Rmax = 5
Nmesh = 6
m = linmesh(Rmin,Rmax,Nmesh)
normalize = true
ordermin = 2
ordermax = 3
left = false
right = false

basis = ShortP1IntLegendreBasis(m, T; ordermin = ordermin, ordermax = ordermax, left = left, right = right, normalize = normalize)

# Plots of basis
X = LinRange(Rmin, Rmax, Nmesh * 100)
plt_basis = plot()
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(plt_basis, X, p.(X))
end
plt_basis



# Plots of the derivatives of the basis
deriv_basis = deriv(basis)
X = LinRange(Rmin, Rmax, Nmesh * 100)
plt_derivbasis = plot()
for i ∈ 1:length(deriv_basis)
    p = build_basis(deriv_basis, i)
    plot!(plt_derivbasis, X, p.(X))
end
plt_derivbasis


# Assembly matrix
M₀  = mass_matrix(basis)
M₋₁ = weight_mass_matrix(basis, -1)
#M₋₂ = weight_mass_matrix(basis, -2)
A   = mass_matrix(deriv(basis))


