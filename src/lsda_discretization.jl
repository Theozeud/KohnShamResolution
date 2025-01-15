#####################################################################
#                          LSDA Cache
#####################################################################

mutable struct LSDACache
    A               # Matrix of Qᵢ'Qⱼ'
    M₀              # Matrix of QᵢQⱼ
    M₋₁             # Matrix of 1/x QᵢQⱼ
    M₋₂             # Matrix of 1/x^2 QᵢQⱼ
    F               # Tensor of 1/x QᵢQⱼQₖ
    B               # Matrix of 4πρQᵢ
    C               # Matrix for V(ρ) - solution of the Gauss Electrostatic law
    Cᵨ              # Number equal to total charge
    Kin             # Matrix VF of Kinetic 
    Coulomb         # Matrix VF of Coulomb
    Hfix            # Part of Hamilotnian not needing to be recomputed (Kinetic + Colombial)         
    Hartree         # Matrix VF of hartree 
    Vxc             # Matrix VF of UP   Echange-correlation potential
    Energy          # Total Energy
    Energy_kin      # Kinetic Energy
end

mutable struct LSDA_tmp_Cache
    tmp_H           # Store Hₗ
    tmp_D           # Store ρ
    tmp_Dstar       # Store ρ
    tmp_U           # Store the eigenvector
    tmp_ϵ           # Store the eigenvalue
    tmp_n           # Store the occupation number
    tmp_MV          # Store the contraction F:C  
    tmp_index_sort
end

function create_cache_lsda(lₕ, Nₕ, T)
    A           = zeros(T, Nₕ, Nₕ) 
    M₀          = zeros(T, Nₕ, Nₕ)
    M₋₁         = zeros(T, Nₕ, Nₕ)
    M₋₂         = zeros(T, Nₕ, Nₕ)
    F           = zeros(T, Nₕ, Nₕ, Nₕ)
    B           = zeros(T, Nₕ)
    C           = zeros(T, Nₕ)
    Cᵨ          = zero(T)
    Kin         = zeros(T, lₕ+1, Nₕ, Nₕ)
    Coulomb     = zeros(T, lₕ+1, Nₕ, Nₕ)
    Hfix        = zeros(T, lₕ+1, Nₕ, Nₕ)
    Hartree     = zeros(T, Nₕ, Nₕ)
    Vxc         = zeros(T, 2, Nₕ, Nₕ)
    Energy      = zero(T)
    Energy_kin  = zero(T)

    # Initialization of array for temporary stockage of computations
    tmp_H           = zeros(T, 2, lₕ+1, Nₕ, Nₕ)
    tmp_D           = zeros(T, 2, Nₕ, Nₕ)
    tmp_Dstar       = zeros(T, 2, Nₕ, Nₕ)
    tmp_U           = zeros(T, 2, lₕ+1, Nₕ, Nₕ)
    tmp_MV          = zeros(T, Nₕ, Nₕ)
    tmp_ϵ           = zeros(T, 2, lₕ+1, Nₕ)
    tmp_n           = zeros(T, 2, lₕ+1, Nₕ)   
    tmp_index_sort  = zeros(Int, Nₕ*(lₕ+1))
    LSDACache(A, M₀, M₋₁, M₋₂, F, B, C, Cᵨ, Kin, Coulomb, Hfix, Hartree, Vxc, Energy, Energy_kin),  
    LSDA_tmp_Cache(tmp_H, tmp_D, tmp_Dstar, tmp_U,  tmp_ϵ, tmp_n, tmp_MV, tmp_index_sort)
end


#####################################################################
#                          LSDA Discretization
#####################################################################

struct LSDADiscretization{T} <: KohnShamDiscretization
    lₕ::Int
    Nₕ::Int
    basis::Basis
    mesh::Mesh{T}
    Rmin::T
    Rmax::T
    elT::Type
    cache::LSDACache
    tmp_cache::LSDA_tmp_Cache
    function LSDADiscretization(lₕ::Int, basis::Basis, mesh::Mesh)
        elT = try
             bottom_type(basis)
        catch
            eltype(basis)
        end
        Nₕ = length(basis)
        new{eltype(mesh)}(lₕ, Nₕ, basis, mesh, first(mesh), last(mesh), elT, create_cache_lsda(lₕ, Nₕ, elT)...)
    end
end

#####################################################################
#                          Init Cache
#####################################################################

function init_cache!(discretization::LSDADiscretization, model::AbstractDFTModel, hartree)

    @unpack lₕ, basis  = discretization
    @unpack A, M₀, M₋₁, M₋₂, F, Kin, Coulomb, Hfix = discretization.cache

    # Creation of the base matrices
    A   .= stiffness_matrix(basis)
    M₀  .= mass_matrix(basis)
    M₋₁ .= weight_mass_matrix(basis, -1)

    # Creation of the fix part of the hamiltonian   
    if lₕ ≠ 0
        M₋₂ .= weight_mass_matrix(basis, -2)
    end
    kinetic_matrix!(discretization)
    coulomb_matrix!(discretization, model)
    @. Hfix = Kin + Coulomb

    # Creation of the 3-index tensor F if there is the hartree term
    if !iszero(hartree)
        F .= weight_mass_3tensor(basis, Monomial(-1))
    end

    nothing
end


#####################################################################
#                          Init for Solver
#####################################################################

init_density_matrix(kd::LSDADiscretization)        = zeros(kd.elT, 2, kd.Nₕ, kd.Nₕ)  
init_coeffs_discretization(kd::LSDADiscretization) = zeros(kd.elT, 2, kd.lₕ+1, kd.Nₕ, kd.Nₕ)
init_energy(kd::LSDADiscretization)                = zeros(kd.elT, 2, kd.lₕ+1, kd.Nₕ)
init_occupation(kd::LSDADiscretization)            = zeros(kd.elT, 2, kd.lₕ+1, kd.Nₕ)

#####################################################################
#               Find Orbital : Solve the eigen problems
#####################################################################

function find_orbital!(discretization::LSDADiscretization, solver::KhonShamSolver)

    @unpack lₕ = discretization
    @unpack M₀, Hfix, Hartree, Vxc = discretization.cache
    @unpack tmp_H, tmp_U, tmp_ϵ = discretization.tmp_cache
    @unpack Dprev = solver
    @unpack hartree = solver.opts

    # STEP 1 : Compute Hartree term 
    if !iszero(hartree)
        hartree_matrix!(discretization, Dprev)
        @. Hartree = hartree * Hartree
    end

    # STEP 2 : Compute Exchange Correlation term
    if isthereExchangeCorrelation(solver.model)
        exchange_corr_matrix!(discretization, solver.model, Dprev)
    end

    # STEP 3 : Solve the generalized eigenvalue problem for each section l and spin up and down
    @threads for σ ∈ 1:2
        @threads for l ∈ 0:lₕ
            # building the hamiltonian of the lᵗʰ section
            @. tmp_H[σ, l+1,:,:] = Hfix[l+1,:,:] + Vxc[σ,:,:] + Hartree
            # solving
            tmp_ϵ[σ, l+1,:], tmp_U[σ, l+1,:,:] = solve_generalized_eigenvalue_problem(tmp_H[σ, l+1,:,:], M₀)
        end
    end
end

#####################################################################
#               Normamisation of eigenvector
#####################################################################

function normalization!(discretization::LSDADiscretization)
    @unpack M₀ = discretization.cache
    @unpack tmp_U, tmp_n = discretization.tmp_cache
    @unpack lₕ, Nₕ  = discretization
    @inbounds for σ ∈ 1:2
        @inbounds for l ∈ 1:lₕ+1
            @inbounds for k ∈ 1:Nₕ
                if !iszero(tmp_n[l,k])
                    normalization = sqrt(tmp_U[σ,l,:,k]' * M₀ * tmp_U[σ,l,:,k])
                    tmp_U[σ,l,:,k] .= tmp_U[σ,l,:,k] .* 1/normalization
                end
            end
        end
    end
    nothing
end

#####################################################################
#                          Kinetic Matrix
#####################################################################

function kinetic_matrix!(discretization::LSDADiscretization)
    @unpack A, M₋₂, Kin = discretization.cache
    for l ∈ 0:discretization.lₕ
        @. Kin[l+1,:,:] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end

#####################################################################
#                          Coulomb Matrix
#####################################################################

function coulomb_matrix!(discretization::LSDADiscretization, model)
    @unpack M₋₁, Coulomb = discretization.cache
    for l ∈ 0:discretization.lₕ
        Coulomb[l+1,:,:] .= - model.z .* M₋₁
    end 
    nothing
end

#####################################################################
#                          Hartree Matrix
#####################################################################

function hartree_matrix!(discretization::LSDADiscretization, D)
    @unpack A, M₀, F, B, C, Hartree = discretization.cache
    @unpack tmp_MV = discretization.tmp_cache
    @unpack basis, Rmin, Rmax = discretization
    @views DUP = D[1,:,:]   
    @views DDOWN = D[2,:,:]
    DD = DUP .+ DDOWN 
    @tensor B[m] = DD[i,j] * F[i,j,m]
    C .= A\B
    @tensor newCᵨ = DD[i,j] * M₀[i,j]
    @tensor tmp_MV[i,j] = C[k] * F[i,j,k]
    @. Hartree = tmp_MV + newCᵨ/(Rmax-Rmin) * M₀
    discretization.cache.Cᵨ = newCᵨ
    nothing
end

#####################################################################
#                   Exchange Correlation Matrix
#####################################################################

function exchange_corr_matrix!(discretization::LSDADiscretization, model, D)
    @views DUP = D[1,:,:]   
    @views DDOWN = D[2,:,:]
    exchange_corr_matrixUP!(discretization, model, DUP)
    exchange_corr_matrixDOWN!(discretization, model, DDOWN)
    nothing
end

function exchange_corr_matrixUP!(discretization::LSDADiscretization, model, DUP)
    @unpack Vxc = discretization.cache
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityUP(discretization, DDOWN, x)
    weightUP(x) = vxcUP(model.exc, ρUP(x), ρDOWN(x))
    fill_weight_mass_matrix!(Vxc[1,:,:], discretization.basis, weightUP)
    nothing
end

function exchange_corr_matrixDOWN!(discretization::LSDADiscretization, model, DDOWN)
    @unpack Vxc = discretization.cache
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityUP(discretization, DDOWN, x)
    weightDOWN(x) = vxcDOWN(model.exc, ρUP(x), ρDOWN(x))
    fill_weight_mass_matrix!(Vxc[2,:,:], discretization.basis, weightDOWN)
    nothing
end

#####################################################################
#                             Energy
#####################################################################

function compute_energy!(discretization::LSDADiscretization)
    compute_total_energy!(discretization)
    compute_kinetic_energy!(discretization)
end

function compute_total_energy!(discretization::LSDADiscretization)
    @unpack Rmin, Rmax = discretization
    @unpack B, C, Cᵨ = discretization.cache
    @unpack tmp_n, tmp_ϵ = discretization.tmp_cache
    @tensor energy = tmp_n[σ, l,n] * tmp_ϵ[σ, l,n] 
    discretization.cache.Energy = energy - discretization.elT(0.5) * (dot(B,C) + Cᵨ^2/(Rmax-Rmin))
    nothing
end

function compute_kinetic_energy!(discretization::LSDADiscretization)

    nothing
end


#####################################################################
#                             Density
#####################################################################

function density_matrix!(discretization::LSDADiscretization)
    @unpack tmp_Dstar, tmp_U, tmp_n = discretization.tmp_cache
    @unpack lₕ, Nₕ  = discretization
    @inbounds for σ ∈ 1:2
        @inbounds for l ∈ 1:lₕ+1 
            @inbounds for k ∈ 1 :Nₕ
                if !iszero(tmp_n[σ, l,k])
                    @inbounds for i ∈ 1:Nₕ
                        val = tmp_n[σ,l,k] * tmp_U[σ,l,i,k] 
                        @inbounds @simd for j ∈ 1:i
                            tmp_Dstar[σ,i,j] += val * tmp_U[σ,l,j,k]
                        end
                    end
                end
            end
        end
    end
    @inbounds for i in 1:Nₕ
        @inbounds @simd for j in 1:i-1
            tmp_Dstar[:,j,i] = tmp_Dstar[:,i,j]
        end
    end
    nothing
end

function compute_density(discretization::LSDADiscretization, D, x)
    @unpack basis = discretization
    newT = promote_type(eltype(basis), typeof(x))
    val = zero(newT)
    eval_basis = zeros(newT, basis.size)
    @inbounds for i ∈ eachindex(basis)
        eval_basis[i] = basis(i,x)
    end
    @views DUP = D[1,:,:]   
    @views DDOWN = D[2,:,:]
    val = (eval_basis)' * (DUP + DDOWN) * eval_basis
    return val* 1/4π * 1/(x^2)
end

function compute_densityUP(discretization::LSDADiscretization, DUP, x)
    @unpack basis = discretization
    newT = promote_type(eltype(basis), typeof(x))
    val = zero(newT)
    eval_basis = zeros(newT, basis.size)
    @inbounds for i ∈ eachindex(basis)
        eval_basis[i] = basis(i,x)
    end
    val = (eval_basis)' * DUP * eval_basis
    return val* 1/4π * 1/(x^2)

end

function compute_densityDOWN(discretization::LSDADiscretization, DDOWN, x)
    @unpack basis = discretization
    newT = promote_type(eltype(basis), typeof(x))
    val = zero(newT)
    eval_basis = zeros(newT, basis.size)
    @inbounds for i ∈ eachindex(basis)
        eval_basis[i] = basis(i,x)
    end
    val = (eval_basis)' * DDOWN * eval_basis
    return val* 1/4π * 1/(x^2)
end