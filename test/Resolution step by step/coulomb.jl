using KohnShamResolution
using Plots
using LinearAlgebra

# General Discretization Parameters
T = Float64
Rmin = 0
Rmax = 10
Nmesh = 10
m = linmesh(Rmin,Rmax,Nmesh)
normalize = true

# Theoretical Eigenvalue and Eigenvectors
f = -Polynomial([Rmin*Rmax, -Rmin - Rmax, 1], 2)
f_true = Monomial(1) * f

# With P1
left = false
right = false
basis = ShortP1Basis(m, T; left = left, right = right, normalize = normalize)
CP1   = weight_mass_matrix(basis, -1)
F = [integrate(f * build_basis(basis, i), Rmin, Rmax) for i ∈ eachindex(basis)]
sol_p1 = build_on_basis(basis, CP1\F)


# With P1-Integrated Legendre Polynomials ordre 2
ordermin = 2
ordermax = 2
basis = ShortP1IntLegendreBasis(m, T; ordermin = ordermin, ordermax = ordermax,  normalize = normalize, left = false, right = false)
CIL2   = weight_mass_matrix(basis, -1)
F = [integrate(f * build_basis(basis, i), Rmin, Rmax) for i ∈ eachindex(basis)]
sol_il2 = build_on_basis(basis, CIL2\F)


# With P1-Integrated Legendre Polynomials ordre 3
ordermin = 2
ordermax = 3
basis = ShortP1IntLegendreBasis(m, T; ordermin = ordermin, ordermax = ordermax,  normalize = normalize)
CIL3   = weight_mass_matrix(basis, -1)
F = [integrate(f*build_basis(basis, i), Rmin, Rmax) for i ∈ eachindex(basis)]
sol_il3 = build_on_basis(basis, CIL3\F)


# Plot solution
X = LinRange(Rmin, Rmax, Nmesh * 100)
plt = plot()

plot!(plt, X, sol_p1.(X), label = "p1", lw = 3)
plot!(plt, X, sol_il2.(X), label = "Integrated Legendre ordre 2", lw = 3)
plot!(plt, X, sol_il3.(X), label = "Integrated Legendre ordre 3", lw = 3)
plot!(plt, X, f_true.(X), label = "Theoretical", lw = 2, color = :black, ls = :dash)

plt