using KohnShamResolution
using Plots

# Test on Bubble basis
T = Float64
cutting_pre = 10
Rmin = 0
Rmax = cutting_pre*log(10)
m = logmesh(Rmin, Rmax, 150)

pb = BubbleBasis(m, T; order = 3, left = false, right = false)


X = LinRange(Rmin,Rmax,1000)
plt = plot(legend = false)
for p ∈ pb
    plot!(X,p.(X))
end
plt

#A = mass_matrix(pb,1,5)