using KohnShamResolution
using Plots
# Choice of the method
method = ConstantODA(1.0)

# Discretization 
lₕ = 0
Rmin = 0.0000001
cutting_pre = 10
T = Float64
ordermin = 2
ordermax = 4

# One electron model
z = 1
N = 1
eigvalue_theo(n,z) = -z^2/(2*n^2)

# Plot
Nmax = 4
Nerror = []


cutting_pre = 30

for n in 2:Nmax

    Rmax = 100#(1.5 * log(z) + cutting_pre*log(10))/z
    m = logmesh(Rmin, Rmax, 2^n; z = 1)
    p1 = ShortP1Basis(m, T;  normalize = true, left = false, right = false)
    intleg = ShortIntLegendreBasis(m, T; normalize = true, ordermin = ordermin, ordermax = ordermax)
    basis = CombineShortPolynomialBasis(p1, intleg)
    D = KohnShamSphericalDiscretization(lₕ, basis, m)
    KM = KohnShamExtended(z = z, N = N)

    @time sol = groundstate(KM, D, method; tol = 1e-20, hartree = false)

    index_ϵ = findall(x->x < 0, sol.ϵ)

    true_ϵ = eigvalue_theo.(index_ϵ, z)

    push!(Nerror, abs.(sol.ϵ[index_ϵ] .- true_ϵ))
end

plterror = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :topright, xaxis=:log, yaxis=:log,
            legendfontsize  = 12,  
            titlefontsize   = 12,
            guidefontsize   = 12,
            tickfontsize    = 12)

xlabel!("Number of point")
ylabel!("Error")

logerror =  []
for i ∈ eachindex(Nerror[end])
    logerror_i = []
    for j ∈ eachindex(Nerror)
        if i ∈ eachindex(Nerror[j])
            push!(logerror_i, Nerror[j][i])
        end
    end
    push!(logerror, logerror_i)
end

for i ∈ eachindex(logerror)
    plot!(2 .^(1+Nmax-length(logerror[i]):Nmax), logerror[i], lw = 3, label = "ϵ"*string(i), 
            markershape = :x, markersize = 10)
end

savefig(plterror, "image/hydrogenoide/with short IntLeg - P1 elements/Convergence")