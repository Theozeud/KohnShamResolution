using KohnShamResolution
using LinearAlgebra
using Plots

# Fixed Parameters 
Rmin = 0.00001
T = Float64

# Parameter of the basis
ordermax = 10
BasisName = "DiffLeg "*string(ordermax)

# Default Parameters
d_Rmax = 100
d_Nmesh = 20
d_z = 1

# Matrix that we check the conditionning
A   = (basis, z) -> mass_matrix(deriv(basis))
M₀  = (basis, z) -> mass_matrix(basis)
M₋₁ = (basis, z) -> weight_mass_matrix(basis, -1)
H   = (basis, z) -> 1/2 * A(basis,z) - z .* M₋₁(basis,z)
fun_Matrix = [A, M₀, M₋₁, H]
label = ["A", "M₀", "M₋₁", "H"]

########################################
# Condition number as a function of Rmax

# Vector of condition number
cn_Rmax = [[], [], [], []]
n_Rmax = [[], [], [], []]

RmaxArray = 20:5:150

for Rmax in RmaxArray

    # Creation of the mesh
    m = logmesh(Rmin, Rmax, d_Nmesh; z = 0.5)

    # Creation of the basis
    p1 = ShortP1Basis(m, T;  normalize = true, left = false, right = false)
    diffleg = ShortDiffLegendreBasis(m, T; normalize = true, ordermax = ordermax)
    basis = CombineShortPolynomialBasis(p1, diffleg)

    # Compute condition number
    for (i, M) ∈ enumerate(fun_Matrix)
        try 
            ϵ, _ = eigen(M(basis, d_z))
            cond = abs(ϵ[end]) > abs(ϵ[begin]) ? abs(ϵ[end]/ϵ[begin]) : abs(ϵ[begin]/ϵ[end])
            push!(cn_Rmax[i], cond)
            push!(n_Rmax[i], Rmax)
        catch
        end
    end
end


# Plots
plt_Rmax = plot(  size = (900,600), margin = 0.5Plots.cm, legend = :topright, 
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)          
xlabel!(plt_Rmax, "Rmax")
ylabel!(plt_Rmax, "Condition number")
title!(plt_Rmax, "z = "*string(d_z)*", Nmesh = "*string(d_Nmesh))

for i ∈ eachindex(fun_Matrix)
    plot!(plt_Rmax, n_Rmax[i], cn_Rmax[i], lw = 3, markershape = :x, markersize = 8, label = label[i], yaxis=:log)
end

########################################
# Condition number as a function of z

# Vector of condition number
cn_z = [[], [], [], []]
n_z = [[], [], [], []]

zArray = 1:1:20

for z in zArray

    # Creation of the mesh
    m = logmesh(Rmin, d_Rmax, d_Nmesh; z = 0.5)

    # Creation of the basis
    p1 = ShortP1Basis(m, T;  normalize = true, left = false, right = false)
    diffleg = ShortDiffLegendreBasis(m, T; normalize = true, ordermax = ordermax)
    basis = CombineShortPolynomialBasis(p1, diffleg)

    # Compute condition number
    for (i, M) ∈ enumerate(fun_Matrix)
        try 
            ϵ, _ = eigen(M(basis, z))
            cond = abs(ϵ[end]) > abs(ϵ[begin]) ? abs(ϵ[end]/ϵ[begin]) : abs(ϵ[begin]/ϵ[end])
            push!(cn_z[i], cond)
            push!(n_z[i], z)
        catch
        end
    end
end


# Plots
plt_z = plot(  size = (900,600), margin = 0.5Plots.cm, legend = :topright, 
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)          
xlabel!(plt_z, "z")
ylabel!(plt_z, "Condition number")
title!(plt_z, "Rmax = "*string(d_Rmax)*", Nmesh = "*string(d_Nmesh))

for i ∈ eachindex(fun_Matrix)
    plot!(plt_z, n_z[i], cn_z[i], lw = 3, markershape = :x, markersize = 8, label = label[i], yaxis=:log)
end

#########################################
# Condition number as a function of Nmesh

# Vector of condition number
cn_Nmesh = [[], [], [], []]
n_Nmesh = [[], [], [], []]

NmeshArray = 10:10:300

for Nmesh in NmeshArray

    # Creation of the mesh
    m = logmesh(Rmin, d_Rmax, Nmesh; z = 0.5)

    # Creation of the basis
    p1 = ShortP1Basis(m, T;  normalize = true, left = false, right = false)
    diffleg = ShortDiffLegendreBasis(m, T; normalize = true, ordermax = ordermax)
    basis = CombineShortPolynomialBasis(p1, diffleg)

    # Compute condition number
    for (i, M) ∈ enumerate(fun_Matrix)
        try 
            ϵ, _ = eigen(M(basis, d_z))
            cond = abs(ϵ[end]) > abs(ϵ[begin]) ? abs(ϵ[end]/ϵ[begin]) : abs(ϵ[begin]/ϵ[end])
            push!(cn_Nmesh[i], cond)
            push!(n_Nmesh[i], Nmesh)
        catch
        end
    end
end


# Plots
plt_Nmesh = plot(  size = (900,600), margin = 0.5Plots.cm, legend = :topright, 
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)          
xlabel!(plt_Nmesh, "Nmesh")
ylabel!(plt_Nmesh, "Condition number")
title!(plt_Nmesh, "Rmax = "*string(d_Rmax)*", z = "*string(d_z))

for i ∈ eachindex(fun_Matrix)
    plot!(plt_Nmesh, n_Nmesh[i], cn_Nmesh[i], lw = 3, markershape = :x, markersize = 8, label = label[i], yaxis=:log)
end

## Plots

pltfin = plot(plt_Rmax, plt_z, plt_Nmesh, layout = (2,2), size = (2400,1000), margin = 1.2Plots.cm)
pltfin

savefig(pltfin, "image/hydrogenoide/with short diff leg/Condition number pour ordremax = "*string(ordermax))