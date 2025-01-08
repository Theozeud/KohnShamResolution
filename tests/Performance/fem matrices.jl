using KohnShamResolution
using LinearAlgebra
using BenchmarkTools

original_stdout = stdout
output_file = open("tests/Performance/fem.txt", "w")
redirect_stdout(output_file)

# IntLeg2 - Nmesh = 50
println("IntLeg2, Nmesh = 50")
Nmesh = 50
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, N; T = T)
basis = ShortP1IntLegendreBasis(m, T; normalize = false, ordermax = 2)
deriv_basis = deriv(basis)

println("A")
@time A   = mass_matrix(deriv_basis)
println("M0")
@time M₀  = mass_matrix(basis)
println("M₋₁")
@time M₋₁ = weight_mass_matrix(basis, -1)
println("M₋₂")
@time FM₋₂ = weight_mass_matrix(basis, -2)


redirect_stdout(original_stdout)
close(output_file)

println("Performance finished")