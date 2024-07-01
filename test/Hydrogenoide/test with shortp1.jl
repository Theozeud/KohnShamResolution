using KohnShamResolution
using LinearAlgebra
using Plots

# Creation of the model
z = 16
N = 1
KM = KohnShamExtended(z = z,N = N)

# Choice of the method
method = ConstantODA(1.0)

# Discretization 
Nmesh = 200
lₕ = 0
Rmin = 0.00000001
cutting_pre = 10
Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
m = logmesh(Rmin, Rmax, Nmesh; z = 0.5)

T = Float64
basis = P1Basis(m, T; left = false, right = false)
normalize = true
shortbasis = ShortP1Basis(m, T;  normalize = normalize, left = false, right = false)

# Solve
deriv_basis = deriv(basis)
A = mass_matrix(deriv_basis)
M₀ = mass_matrix(basis)
M₋₁ = weight_mass_matrix(basis, -1)
M₋₂ = weight_mass_matrix(basis, -2)

H = 1/2 * (A + lₕ*(lₕ+1)*M₋₂) - z .* M₋₁
ϵ, U = eigen(H,M₀)

shortderiv_basis = deriv(shortbasis)
A = mass_matrix(shortderiv_basis)
M₀ = mass_matrix(shortbasis)
M₋₁ = weight_mass_matrix(shortbasis, -1)
M₋₂ = weight_mass_matrix(shortbasis, -2)

H = 1/2 * (A + lₕ*(lₕ+1)*M₋₂) - z .* M₋₁
ϵshort, Ushort = eigen(H,M₀)

X = LinRange(Rmin,Rmax,1000)
plt = plot()
for i ∈ 1:length(deriv_basis)
    p = build_basis(shortderiv_basis, i)
    plot!(X, p.(X))
end
plt

plot(ϵ)
plot!(ϵshort)