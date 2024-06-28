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
Nmesh = 20
lₕ = 0
Rmin = 0.00000001
cutting_pre = 10
@show (1.5 * log(z) + cutting_pre*log(10))/z
Rmax = 100 #
m = linmesh(Rmin, Rmax, Nmesh)

T = Float64
normalize = true

# Block P1
left = false
right = true
p1 = ShortP1Basis(m, T; normalize = normalize, left = left, right = right)

# Block IntLegendre
ordermin = 2
ordermax = 5
intleg = ShortIntLegendreBasis(m, T; normalize = normalize, ordermin = ordermin, ordermax = ordermax)

# Combinaison
basis = CombineShortPolynomialBasis(p1, intleg)

# Solve
deriv_basis = deriv(basis)
A = mass_matrix(deriv_basis)
M₀ = mass_matrix(basis)
M₋₁ = weight_mass_matrix(basis, -1)
M₋₂ = weight_mass_matrix(basis, -2)

H = 1/2 * (A + lₕ*(lₕ+1)*M₋₂) - z .* M₋₁
ϵ, U = eigen(H,M₀)

#=
X = LinRange(Rmin,Rmax,1000)
plt = plot()
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(X, p.(X))
end
plt
=#

scatter(ϵ)