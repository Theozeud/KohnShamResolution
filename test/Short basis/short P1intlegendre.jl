using KohnShamResolution
using Plots

# Parameters of the discretization
using Quadmath
T = Float128
Rmin = 0
Rmax = 36
Nmesh = 8
m = logmesh(Rmin,Rmax,Nmesh)
normalize = true
ordermin = 2
ordermax = 2
left = false
right = false

basis = ShortP1IntLegendreBasis(m, T; ordermin = ordermin, ordermax = ordermax, left = left, right = right, normalize = normalize, Rcut = 37)


# Plots elements of the basis
X = LinRange(-1, 1, 10000)
plt_elements = plot(legend = false)
for b ∈ basis.basisVector
    for i ∈ eachindex(b.elements)
        plot!(plt_elements, X, b.elements[i].(X))
    end
end
plt_elements


# Plots of basis
X = LinRange(Rmin, Rmax, 10000)
plt_basis = plot(legend = false)
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(plt_basis, X, p.(X))
end
plt_basis


# Plots of the derivatives of the basis
deriv_basis = deriv(basis)
X = LinRange(Rmin, Rmax, 10000)
plt_derivbasis = plot(legend = false)
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


