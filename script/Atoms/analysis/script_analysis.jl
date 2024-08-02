using KohnShamResolution
using LinearAlgebra
using Plots
using LambertW
include("../benchmark tools/firstatoms.jl")

# Parameters
T = Float64
z = 10
N = 10
Rmin = 0
Rmax = Int(round(-lambertw(-1e-16/sqrt(z), -1)/z))
Nmesh = 70
lₕ = 1
maxiter = 10
oda = 0.7
tol = 1e-7
sols, title, label = compute_sol(;z = z, N = N, Rmax = Rmax, Nmesh = Nmesh, lₕ = lₕ, maxiter = maxiter, oda = oda, tol = tol, T = T)

plt_crit = plot_crit(sols, tol; title, label)


plt_Ehisto, plt_diffE = plot_Ehisto(sols; title, label)

plt_density = plot_density(sols, Rmax; title, label)

integrate(sols[3].ρ*Monomial(2), Rmin, Rmax)*4π