using KohnShamResolution
using LinearAlgebra
using BenchmarkTools

original_stdout = stdout
output_file = open("tests/Performance/fem_$(Threads.nthreads())_threads.txt", "w")
redirect_stdout(output_file)



# IntLeg2 - Nmesh = 50
println("IntLeg2, Nmesh = 50")
Nmesh = 50
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, Nmesh; T = T)
println("basis")
@time basis = ShortIntLegendreBasis(m, T; ordermax = 2)

println("A")
@time A   = stiffness_matrix(basis)
println("M0")
@time M₀  = mass_matrix(basis)
println("M₋₁")
@time M₋₁ = weight_mass_matrix(basis, -1)
println("M₋₂")
@time FM₋₂ = weight_mass_matrix(basis, -2)

println("=====================================")

# IntLeg5 - Nmesh = 50
println("IntLeg5, Nmesh = 50")
Nmesh = 50
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, Nmesh; T = T)
println("basis")
@time basis = ShortIntLegendreBasis(m, T; ordermax = 5)

println("A")
@time A   = stiffness_matrix(basis)
println("M0")
@time M₀  = mass_matrix(basis)
println("M₋₁")
@time M₋₁ = weight_mass_matrix(basis, -1)
println("M₋₂")
@time FM₋₂ = weight_mass_matrix(basis, -2)

println("=====================================")

# IntLeg5 - Nmesh = 300
println("IntLeg5, Nmesh = 300")
Nmesh = 300
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, Nmesh; T = T)
println("basis")
@time basis = ShortIntLegendreBasis(m, T; ordermax = 5)

println("A")
@time A   = stiffness_matrix(basis)
println("M0")
@time M₀  = mass_matrix(basis)
println("M₋₁")
@time M₋₁ = weight_mass_matrix(basis, -1)
println("M₋₂")
@time FM₋₂ = weight_mass_matrix(basis, -2)


redirect_stdout(original_stdout)
close(output_file)

println("Performance finished")