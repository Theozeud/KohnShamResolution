using KohnShamResolution

include("../plot tools.jl")

# Choice of the method
method = ODA()

# Discretization 
lₕ = 0
m = logmesh(0,50,100)
basis = HatBasis(m; left = false, right = false)
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# One electron model
zA = [1, 2, 3, 4]
N = 1
eigvalue_theo(n,z) = -z^2/(2*n^2)
fundamental(z, x)     = z^(3/2)/sqrt(π)*exp(-z*abs(x))

# Plot
pltA = []

for z in zA

    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :bottomright,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)

    xlabel!("n")
    ylabel!("Energie")
    title!("z = "*string(z))

    KM = KohnShamExtended(z = z,N = N)

    sol = groundstate(KM, D, method; tol = 1e-20, hartree = false)

    index_ϵ = findall(x->x < 0, sol.ϵ[1,:])

    scatter!(index_ϵ, eigvalue_theo.(index_ϵ,z),  label = "Théorique",
                markershape = :circ, 
                markersize = 8,
                markeralpha = nothing,
                markercolor = :red,
                markerstrokewidth = 1,
                markerstrokecolor = :red
                )

    scatter!(index_ϵ, sol.ϵ[index_ϵ],  label = "Numérique",
                markershape = :x, 
                markersize = 10,
                markeralpha = nothing,
                markercolor = :black,
                markerstrokewidth = 2)

    push!(pltA, plt)

end

pltfin = plot(pltA..., layout = (2,2), size = (1200,1000))
savefig(pltfin, "Comparaison Numérique - Théorique pour modèle Hartree Fock à un électron.")