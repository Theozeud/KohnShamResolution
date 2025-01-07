using KohnShamResolution
using Test
using Plots


T = Float64

for n ∈ 2:100
    for m ∈ 2:n-1
        P = intLegendre(n)
        Q = intLegendre(m)
        if !(scalar_product(P,Q,-1,1) ≈ 0)
            println((n,m))
        end
    end
end
