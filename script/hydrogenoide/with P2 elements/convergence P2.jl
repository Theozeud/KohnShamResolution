using KohnShamResolution
using Plots

# Choice of the method
method = ODA()

# Discretization 
lₕ = 0

# One electron model
z = 1
N = 1
eigvalue_theo(n,z) = -z^2/(2*n^2)

# Plot
Nmax = 10
Nerror = []

cutting_pre = 10

for n in 2:Nmax

    Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
    m = logmesh(0, Rmax, 2^n; z = 1/z)
    basis = P2Basis(m; left = false, right = false)
    D = KohnShamSphericalDiscretization(lₕ, basis, m)

    KM = KohnShamExtended(z = z, N = N)

    @time sol = groundstate(KM, D, method; tol = 1e-20, hartree = false)

    index_ϵ = findall(x->x < 0, sol.ϵ)

    true_ϵ = eigvalue_theo.(index_ϵ, z)

    push!(Nerror, abs.(sol.ϵ[index_ϵ] .- true_ϵ))
end

plterror = plot(size = (1000,800), margin = 0.5Plots.cm, legend = :outertopright, xaxis=:log, yaxis=:log,
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
    @show logerror[i]
    @show length(logerror[i])
    @show length((1+Nmax-length(logerror[i]):Nmax))
    plot!(2 .^(1+Nmax-length(logerror[i]):Nmax), logerror[i], lw = 3, label = "ϵ"*string(i), 
            markershape = :x, markersize = 10)
end

savefig(plterror, "image/hydrogenoide/with P2 elements/Convergence P2")