using KohnShamResolution
using Plots

# One electron model
zA = [1, 4, 8, 16]
N = 1
eigvalue_theo(n,z) = -z^2/(2*n^2)

# Choice of the method
method = ConstantODA(1.0)

# Discretization 
lₕ = 0
Rmin = 0
Rmax = 100
cutting_pre = 10
Nmesh = 100
T = Float64
ordermin = 2
ordermax = 5

# Plots
pltA = []

for z in zA

    #Rmax = (1.5 * log(z) + cutting_pre*log(10))/z
    m = logmesh(Rmin, Rmax, Nmesh; z = 0.5)

    p1 = ShortP1Basis(m, T;  normalize = true, left = false, right = false)
    intleg = ShortIntLegendreBasis(m, T; normalize = true, ordermin = ordermin, ordermax = ordermax)
    basis = CombineShortPolynomialBasis(p1, intleg)

    D = KohnShamSphericalDiscretization(lₕ, basis, m)
    KM = KohnShamExtended(z = z, N = N)

    sol = groundstate(KM, D, method; tol = 1e-4, hartree = false)

    plt = plot( size = (900,600), margin = 0.5Plots.cm, legend = :bottomright,
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)

    xlabel!("n")
    ylabel!("Energie")
    title!("z = "*string(z))

    index_ϵ = findall(x->x < 0, sol.ϵ)

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
savefig(pltfin, "image/hydrogenoide/with short IntLeg - P1 elements/Valeurs propres avec Nmesh = "*string(Nmesh)*" et ordres = "*string(ordermax)*" et Rmax = "*string(Rmax))