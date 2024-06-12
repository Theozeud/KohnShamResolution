using KohnShamResolution
using Plots

# Choice of the method
method = ODA()

# Discretization 
lₕ = 0

function lin_to_log(X)
    Y = sort(X)
    a = first(Y)
    b = last(Y)
    Z = zero(X)
    for i ∈ firstindex(Y):lastindex(Y)-1
        x = Y[i]
        tmp = (b-a)/(b-x)
        tmp = log(tmp)
        Z[i] = b - (b-a)/(tmp+1)
    end
    Z[end] = b
    Z
end

m = mesh(lin_to_log(LinRange(0,50,100)))
basis = HatBasis(m; left = false, right = false)
D = KohnShamSphericalDiscretization(lₕ, basis, m)

# Model(s)

zA = [0.2,0.5,1.0,1.5]
N = 1
theorique(n,z) = -z^2/(2*n^2)

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

    scatter!(index_ϵ, theorique.(index_ϵ,z),  label = "Théorique",
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
