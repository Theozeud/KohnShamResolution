using KohnShamResolution
using LinearAlgebra
using Plots

Rmin = 0
Rmax = 20
N = 20
s = 0.9



function geometricrange2(a,b,n; T = Float64, s)
    R = zeros(T,n)
    R[1] = a
    hn = (one(T)-T(s))/(one(T) - T(s)^n)*(b-a)
    H = zeros(T,n-1)
    @show H[end] = hn
    for i ∈ n-1:-1:2
        @show H[i-1] = T(s) * H[i]
    end 
    @show H
    for i ∈ 2:n
        R[i] = R[i-1] + H[i-1]
    end
    R
end

m = geometricrange2(Rmin, Rmax, N; s = s)

plot(m, markershape = :x)