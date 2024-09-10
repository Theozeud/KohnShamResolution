using KohnShamResolution
using Plots

# Parameters of the discretization
T = Float64
Rmin = 0
Rmax = 5
Nmesh = 90
m = linmesh(Rmin,Rmax,Nmesh)
normalize = true
ordermax = 5

basis = ShortP1IntLegendreBasis(m, T; ordermax = ordermax, normalize = normalize)
deriv_basis = deriv(basis)

# Plots elements of the basis
X = LinRange(-1, 1, 10000)
plt_elements = plot(legend = false)
for b ∈ basis.basisVector
    for i ∈ eachindex(b.elements)
        plot!(plt_elements, X, b.elements[i].(X), lw = 4)
    end
end
plt_elements

################
# Plots of basis
X = LinRange(Rmin, Rmax, 10000)
plt_basis = plot(legend = false)
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(plt_basis, X, p.(X), lw = 4)
end
xlabel!("Maillage")
plt_basis

#######################################
# Plots of the derivatives of the basis
X = LinRange(Rmin, Rmax, 10000)
plt_derivbasis = plot(legend = false)
for i ∈ 1:length(deriv_basis)
    p = build_basis(deriv_basis, i)
    plot!(plt_derivbasis, X, p.(X))
end
plt_derivbasis


######################
# Test Assembly matrix
M₀  = mass_matrix(basis)
A   = mass_matrix(deriv_basis)
M₋₁ = weight_mass_matrix(basis, -1)
M₋₂ = weight_mass_matrix(basis, -2)
F   = weight_mass_3tensor(basis, Monomial(-1))




#######################
# Debuggage

function isthereNaN(A)
    for e ∈ A
        isnan(e)
        return true
    end
    false
end

function debuggage_weight_mass_matrix(cb, weight)
    T = bottom_type(cb)

    # First step : test each block of the matrix
    A = zeros(T, cb.size, cb.size)
    index_blocks_problem = CartesianIndex{2}[]
    blocks_problem = []
    for b ∈ KohnShamResolution.getblocks(cb)
        @views ABlock = A[KohnShamResolution.getrangerow(b), KohnShamResolution.getrangecolumn(b)]
        if KohnShamResolution.isdiagonal(b)
            try 
                KohnShamResolution.fill_weight_mass_matrix!(KohnShamResolution.getbasis(cb, KohnShamResolution._getindex(b,1)), weight, ABlock)
                if isthereNaN(ABlock)
                    push!(index_blocks_problem, b.index)
                    push!(blocks_problem, b)
                end
            catch
                push!(index_blocks_problem, b.index)
                push!(blocks_problem, b)
            end
        else
            try
                KohnShamResolution.fill_weight_mass_matrix!(KohnShamResolution.getbasis(cb, KohnShamResolution._getindex(b,1)), KohnShamResolution.getbasis(cb, KohnShamResolution._getindex(b,2)), weight, b.interaction_index, ABlock)
                if isthereNaN(ABlock)
                    push!(index_blocks_problem, b.index)
                    push!(blocks_problem, b)
                end
            catch
                push!(index_blocks_problem, b.index)
                push!(blocks_problem, b)
            end
        end
    end
    println("There is problem with following blocks : $index_blocks_problem")

    # Test inside each blocsk
    info_problem = []
    for b ∈ blocks_problem
        println("Debugg in the $(b.index) block.")
        spb1 = KohnShamResolution.getbasis(cb, KohnShamResolution._getindex(b,1))
        spb2 = KohnShamResolution.getbasis(cb, KohnShamResolution._getindex(b,2))
        for I ∈ b.interaction_index
            for (i,j) ∈ KohnShamResolution.intersection_with_indices(KohnShamResolution.getsegments(spb1, I[1]), KohnShamResolution.getsegments(spb2, I[2]))
                P = KohnShamResolution.getpolynomial(spb1, I[1], i)
                Q = KohnShamResolution.getpolynomial(spb2, I[2], j)
                invϕ = KohnShamResolution.getinvshift(spb1, I[1], i)
                dinvϕ = invϕ[1]
                try
                    wsp = KohnShamResolution.weight_scalar_product(P, Q, weight, spb1.elements.binf, spb1.elements.bsup, invϕ)
                    if isnan(wsp) || isinf(wsp)
                        push!(info_problem, (b.index,I,P,Q, weight, invϕ))
                    end
                catch
                    push!(info_problem, (b.index,I,P,Q, weight, invϕ))
                end
            end
        end
    end
    info_problem
end

info_problem = debuggage_weight_mass_matrix(basis, Monomial(-2))


@show (index,I,P,Q, weight, invϕ) = info_problem[1]
KohnShamResolution.weight_scalar_product(P, Q, weight, -1, 1, invϕ)


#=
nb1 = 0
for I∈CartesianIndices(M₋₁)
    if isnan(M₋₁[I]) || abs(M₋₁[I])==Inf
        global nb1 += 1
    end
end
println("Thre is $nb1 mistake coefficiens of M₋₁.")

nb2 = 0
for I∈CartesianIndices(M₋₂)
    if isnan(M₋₂[I]) || abs(M₋₂[I])==Inf
        nb2 += 1
    end
end
println("Thre is $nb2 mistake coefficiens of M₋₂.")
=#

