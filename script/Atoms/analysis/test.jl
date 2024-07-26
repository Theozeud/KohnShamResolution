using KohnShamResolution
include("../benchmark tools/firstatoms.jl")
# Parameters
T = Float64
z = 1
N = 1
Rmax = 60
Nmesh = 100
lₕ = 0
maxiter = 1
oda = 0.8
tol = 1e-4

sols, plt_criteria = testvalue(;z = z, N = N, Rmax = Rmax, Nmesh = Nmesh, lₕ = lₕ, maxiter = maxiter, oda = oda, tol = tol, T = T) 

for sol ∈ sols
    @show sol.success
    @show sol.E
end