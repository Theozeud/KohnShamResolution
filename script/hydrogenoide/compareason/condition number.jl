using KohnShamResolution
using LinearAlgebra
using Plots

# Choice of the method
method = ODA()

# Discretization 
lₕ = 0
Rmin = 0

# One electron model
z = 1
N = 1
KM = KohnShamExtended(z = z, N = N) 
eigvalue_theo(n,z)    = -z^2/(2*n^2)
fundamental(z, x)     = z^(3/2)/sqrt(π)*exp(-z*abs(x))

# Vector of condition number

condition_number_Array = []
RmaxArray = 10:5:150

# Basis to test
Basis = [HatBasis, P2Basis]
BasisName = ["P1", "P2"]

for basis_fun in Basis
    condition_number = []
    for Rmax in RmaxArray
        m = logmesh(0,Rmax,100)
        basis = basis_fun(m; left = false, right = false)
        D = KohnShamSphericalDiscretization(lₕ, basis, m)

        # Computation by hand
        deriv_basis = deriv(basis)

        A   = mass_matrix(deriv_basis, Rmin, Rmax)
        M₀  = mass_matrix(basis, Rmin, Rmax)
        M₋₁ = weight_mass_matrix(basis, -1, Rmin, Rmax)
        M₋₂ = weight_mass_matrix(basis, -2, Rmin, Rmax)

        H = 1/2 * (A + lₕ*(lₕ+1)*M₋₂) - z .* M₋₁

        ϵ, U = eigen(H,M₀)

        push!(condition_number, abs(ϵ[end]/ϵ[begin]))
    end
    push!(condition_number_Array, condition_number)
end


# Plots

plt_cn = plot(  size = (900,600), margin = 0.5Plots.cm, legend = :topright, 
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)
                
xlabel!(plt_cn, "Rmax")
ylabel!(plt_cn, "Conditionning number")
title!(plt_cn, "z = "*string(z))

for (i,condition_number) ∈ zip(eachindex(RmaxArray),condition_number_Array)
    plot!(RmaxArray, condition_number, lw = 3, markershape = :x, markersize = 8, label = BasisName[i])
end


plt_cn_zoom = plot(  size = (900,600), margin = 0.5Plots.cm, legend = :topright, 
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)
                
xlabel!("Rmax")
ylabel!("Conditionning number")
title!("z = "*string(z))

for (i,condition_number) ∈ zip(eachindex(RmaxArray),condition_number_Array)
    plot!(RmaxArray[2*div(length(RmaxArray),3):end], condition_number[2*div(length(RmaxArray),3):end], lw = 3, markershape = :x, markersize = 8, label = BasisName[i])
end

pltfin = plot(plt_cn, plt_cn_zoom, layout = (1,2), size = (2400,1000), margin = 1.2Plots.cm)
savefig(pltfin, "image/hydrogenoide/compareason/Condition Number")