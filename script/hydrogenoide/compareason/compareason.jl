using KohnShamResolution

# Choice of the method
method = ODA()

# Discretization 
lₕ = 0

# One electron model
z = 1
N = 1
eigvalue_theo(n,z) = -z^2/(2*n^2)

# Plot
Nmax = 6
Nerror_P1 = []
Nerror_P2 = []

cutting_pre = 10

for n in 2:Nmax

    Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
    m = logmesh(0,Rmax,2^n, 1/z)

    basis_P1 = HatBasis(m; left = false, right = false)
    DP1 = KohnShamSphericalDiscretization(lₕ, basis_P1, m)
    
    basis_P2 = P2Basis(m; left = false, right = false)
    DP2 = KohnShamSphericalDiscretization(lₕ, basis_P2, m)

    KM = KohnShamExtended(z = z, N = N)

    sol_P1 = groundstate(KM, DP1, method; tol = 1e-20, hartree = false)
    sol_P2 = groundstate(KM, DP2, method; tol = 1e-20, hartree = false)

    index_ϵ_P1 = findall(x->x < 0, sol_P1.ϵ[1,:])
    index_ϵ_P2 = findall(x->x < 0, sol_P2.ϵ[1,:])

    true_ϵ_P1 = eigvalue_theo.(index_ϵ_P1, z)
    true_ϵ_P2 = eigvalue_theo.(index_ϵ_P2, z)

    push!(Nerror_P1, sol_P1.ϵ[index_ϵ_P1] .- true_ϵ_P1)
    push!(Nerror_P2, sol_P2.ϵ[index_ϵ_P2] .- true_ϵ_P2)
end

plterror = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :bottomleft, xaxis=:log, yaxis=:log,
            legendfontsize  = 12,  
            titlefontsize   = 12,
            guidefontsize   = 12,
            tickfontsize    = 12)

xlabel!("Number of point")
ylabel!("Error")

logerror_P1 =  []
for i ∈ eachindex(Nerror_P1[end])
    logerror_i = []
    for j ∈ eachindex(Nerror_P1)
        if i ∈ eachindex(Nerror_P1[j])
            push!(logerror_i, Nerror_P1[j][i])
        end
    end
    push!(logerror_P1, logerror_i)
end

logerror_P2 =  []
for i ∈ eachindex(Nerror_P2[end])
    logerror_i = []
    for j ∈ eachindex(Nerror_P2)
        if i ∈ eachindex(Nerror_P2[j])
            push!(logerror_i, Nerror_P2[j][i])
        end
    end
    push!(logerror_P2, logerror_i)
end

for i ∈ eachindex(logerror_P1)[1:end-1]
    plot!(2 .^(1+Nmax-length(logerror_P1[i]):Nmax), logerror_P1[i], lw = 2, label = "ϵ"*string(i), 
            markershape = :x, markersize = 10, ls = :dash)
end

for i ∈ eachindex(logerror_P2)[1:end-1]
    plot!(2 .^(1+Nmax-length(logerror_P2[i]):Nmax), logerror_P2[i], lw = 2, label = "ϵ"*string(i), 
            markershape = :+, markersize = 12, ls= :solid)
end

savefig(plterror, "image/hydrogenoide/compareason/Convergence P1-P2")