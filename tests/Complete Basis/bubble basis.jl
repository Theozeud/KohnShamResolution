using KohnShamResolution
using Plots

# Test on Bubble basis
T = Rational{BigInt}
cutting_pre = 10
Rmin = 0
Rmax = cutting_pre*log(10)
m = logmesh(Rmin, Rmax, 10; T = T)

pb = BubbleBasis(m, T; order = 3, left = false, right = false)

X = LinRange(Rmin, Rmax, 1000)
plt = plot(legend = false)
for p ∈ pb
    plot!(X,p.(X))
end
plt
