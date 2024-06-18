using KohnShamResolution
using LinearAlgebra
using Plots

# Creation of the model
z = 4
N = 1

KM = KohnShamExtended(z = z,N = N)
#KM = SlaterXα(z, N)

# Choice of the method
method = ODA()

# Discretization 
lₕ = 2
Rmin = 0
cutting_pre = 10
Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
m = logmesh(Rmin, Rmax, 10)
basis = P1Basis(m; left = false, right = false)
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Solve
#@time sol = groundstate(KM, D, method; tol = 1e-10, hartree = false)

deriv_basis = deriv(basis)
 
A   = mass_matrix(deriv_basis, Rmin, Rmax)
M₀  = mass_matrix(basis, Rmin, Rmax)
M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

U = zeros(lₕ+1, 8, 8)
ϵ = zeros(lₕ+1, 8)

for l ∈ 0:lₕ
    H = 1/2 * (A + l*(l+1)*M₋₂) - z .* M₋₁
    _ϵ, _U = eigen(H,M₀)
    ϵ[l+1,:] .= _ϵ
    U[l+1,:,:] .= _U
end

#= Plot Fundamental
plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :bottomright,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)

xlabel!("r")
ylabel!("Fundamental")
title!("z = "*string(z))

X = logmesh(Rmin,Rmax, 1000).points

fundamental(z, x)     = z^(3/2)/sqrt(π)*exp(-z*abs(x))
plot!(X, fundamental.(z,X),  label = "Théorique", color = :red, lw = 3)

eig1 = KohnShamResolution.build_on_basis(basis,U[:,1]) 
fun =   1/sqrt(4π) * eig1 / sqrt(normL2(eig1, m))  * Monomial(-1)

plot(plt, X, fun.(X), label = "Numérique", color = :blue)
=#

