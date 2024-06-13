using KohnShamResolution
using Test
using Plots

# Test on P2 basis
T = Float16
m = mesh(1:5)

p2 = P2Basis(m, T; left = false)

X = LinRange(1,5,1000)
plt = plot()
for p ∈ p2
    plot!(X,p.(X))
end
plt