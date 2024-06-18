using KohnShamResolution

# Choice of the method
method = ODA()

# Discretization 
lₕ = 0

# One electron model
zA = [1, 4, 8, 16]
N = 1
eigvalue_theo(n,z) = -z^2/(2*n^2)

# Plot
pltA = []

cutting_pre = 10
order = 3
T = Rational{BigInt}

for z in zA

    Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
    m = logmesh(0, Rmax, 15; z=  1/z, T = T)
    basis = BubbleBasis(m; order = order, left = false, right = false)
    D = KohnShamSphericalDiscretization(lₕ, basis, m)

    KM = KohnShamExtended(z = z, N = N)

    sol = groundstate(KM, D, method; tol = 1e-20, hartree = false)

    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :bottomright,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)

    xlabel!("n")
    ylabel!("Energie")
    title!("z = "*string(z))

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
savefig(pltfin, "image/hydrogenoide/with Bubble elements/Valeurs propres Ordre "*string(order))