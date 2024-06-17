using KohnShamResolution
using Plots

# Test on Bubble basis
T = Float16
m = mesh(1:5)

pb = BubbleBasis(m, T; order = 1, left = false)

X = LinRange(1,5,1000)
plt = plot()
for p ∈ pb
    plot!(X,p.(X))
end
plt