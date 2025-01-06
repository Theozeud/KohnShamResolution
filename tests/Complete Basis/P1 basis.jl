using KohnShamResolution
using Test
using Plots

# Test on P2 basis
T = Float64
m = mesh(1:5)

p1 = P1Basis(m, T)

X = LinRange(1,5,1000)
plt = plot()
for p ∈ p1
    plot!(X,p.(X))
end
plt
