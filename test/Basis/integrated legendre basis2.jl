using KohnShamResolution
using Plots

# Test on integrated Legendre basis
T = Float64
Rmin = 1
Rmax = 5
Nmesh = 1000
m = linmesh(Rmin, Rmax, Nmesh)

ord = 3
basis = IntLegendreBasis(m, T; order = ord, left = false, right = true)

X = LinRange(-1, 1,1000)
plt = plot(legend = :outertopright)
for p ∈ basis.polynomial
    @show integrate(p*p,-1,1)
    plot!(X,p.(X))
end
plt

@time mass_matrix(basis)

#weight_mass_matrix(basis, 0)