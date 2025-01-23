#####################################################################
#                          LDA Cache
#####################################################################

mutable struct LDACache
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
    Vxc             # Matrix VF of Echange-correlation
end

mutable struct LDA_tmp_Cache
    tmp_H           # Store Hₗ
    tmp_D           # Store ρ
    tmp_Dstar       # Store ρ*
    tmp_U           # Store the eigenvector
    tmp_ϵ           # Store the eigenvalue
    tmp_n           # Store the occupation number
    tmp_MV          # Store the contraction F:C  
    tmp_index_sort
end

function create_cache_lda(lₕ, Nₕ, T, lmin)
    A           = zeros(T, Nₕ, Nₕ) 
    M₀          = zeros(T, Nₕ, Nₕ)
    M₋₁         = zeros(T, Nₕ, Nₕ)
    M₋₂         = zeros(T, Nₕ, Nₕ)
    F           = zeros(T, Nₕ, Nₕ, Nₕ)
    B           = zeros(T, Nₕ)
    C           = zeros(T, Nₕ)
    Cᵨ          = zero(T)
    Kin         = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    Coulomb     = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    Hfix        = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    Hartree     = zeros(T, Nₕ, Nₕ)
    Vxc         = zeros(T, Nₕ, Nₕ)

    # Initialization of array for temporary stockage of computations
    tmp_H           = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    tmp_D           = zeros(T, Nₕ, Nₕ) 
    tmp_Dstar       = zeros(T, Nₕ, Nₕ) 
    tmp_U           = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    tmp_MV          = zeros(T, Nₕ, Nₕ)
    tmp_ϵ           = zeros(T, lₕ+1 - lmin, Nₕ)
    tmp_n           = zeros(T, lₕ+1 - lmin, Nₕ)   
    tmp_index_sort  = zeros(Int, Nₕ*(lₕ+1 - lmin))
    LDACache(A, M₀, M₋₁, M₋₂, F, B, C, Cᵨ, Kin, Coulomb, Hfix, Hartree, Vxc),  
    LDA_tmp_Cache(tmp_H, tmp_D, tmp_Dstar, tmp_U,  tmp_ϵ, tmp_n, tmp_MV, tmp_index_sort)
end


#####################################################################
#                          LDA Discretization
#####################################################################


struct LDADiscretization{T} <: KohnShamDiscretization
    lmin::Int
    lₕ::Int
    Nₕ::Int
    basis::Basis
    mesh::Mesh{T}
    Rmin::T
    Rmax::T
    elT::Type
    cache::LDACache
    tmp_cache::LDA_tmp_Cache
    function LDADiscretization(lₕ::Int, basis::Basis, mesh::Mesh; lmin = 0)
        elT = try
             bottom_type(basis)
        catch
            eltype(basis)
        end
        Nₕ = length(basis)
        @assert lmin ≤ lₕ
        new{eltype(mesh)}(lmin, lₕ, Nₕ, basis, mesh, first(mesh), last(mesh), elT, create_cache_lda(lₕ, Nₕ, elT, lmin)...)
    end
end

#####################################################################
#                          Init Cache
#####################################################################

function init_cache!(discretization::LDADiscretization, model::AbstractDFTModel, hartree)

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

init_density_matrix(kd::LDADiscretization)        = zeros(kd.elT, kd.Nₕ, kd.Nₕ)  
init_coeffs_discretization(kd::LDADiscretization) = zeros(kd.elT, kd.lₕ+1, kd.Nₕ, kd.Nₕ)
init_energy(kd::LDADiscretization)                = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)
init_occupation(kd::LDADiscretization)            = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)

#####################################################################
#               Find Orbital : Solve the eigen problems
#####################################################################

function find_orbital!(discretization::LDADiscretization, solver::KhonShamSolver)

    @unpack lmin, lₕ = discretization
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

    # STEP 3 : Solve the generalized eigenvalue problem for each section l
    @threads for l ∈ lmin:lₕ
        # building the hamiltonian of the lᵗʰ section
        @. tmp_H[l+1-lmin,:,:] = Hfix[l+1-lmin,:,:] + Vxc + Hartree
        # solving
        tmp_ϵ[l+1-lmin,:], tmp_U[l+1-lmin,:,:] = solve_generalized_eigenvalue_problem(tmp_H[l+1-lmin,:,:], M₀)
        # normalization of eigenvector
        
    end
end

#####################################################################
#               Normamisation of eigenvector
#####################################################################

function normalization!(discretization::LDADiscretization)
    @unpack M₀ = discretization.cache
    @unpack tmp_U, tmp_n = discretization.tmp_cache
    @unpack lₕ, Nₕ  = discretization
    @inbounds for l ∈ 1:lₕ+1   
        @inbounds for k ∈ 1:Nₕ
            if !iszero(tmp_n[l,k])
                normalization = sqrt(sum([tmp_U[l,i,k] * tmp_U[l,j,k] * M₀[i,j] for i∈1:Nₕ for j∈1:Nₕ]))
                tmp_U[l,:,k] .= tmp_U[l,:,k] .* 1/normalization
            end
        end
    end
    nothing
end

#####################################################################
#                          Kinetic Matrix
#####################################################################

function kinetic_matrix!(discretization::LDADiscretization)
    @unpack A, M₋₂, Kin = discretization.cache
    for l ∈ discretization.lmin:discretization.lₕ
        @. Kin[l+1-discretization.lmin,:,:] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end

#####################################################################
#                          Coulomb Matrix
#####################################################################

function coulomb_matrix!(discretization::LDADiscretization, model)
    @unpack M₋₁, Coulomb = discretization.cache
    for l ∈ discretization.lmin:discretization.lₕ
        Coulomb[l+1-discretization.lmin,:,:] .= - model.z .* M₋₁
    end 
    nothing
end

#####################################################################
#                          Hartree Matrix
#####################################################################

function hartree_matrix!(discretization::LDADiscretization, D)
    @unpack A, M₀, F, B, C, Hartree = discretization.cache
    @unpack tmp_MV = discretization.tmp_cache
    @unpack basis, Rmin, Rmax = discretization
    @tensor B[m] = D[i,j] * F[i,j,m]
    C .= A\B
    @tensor newCᵨ = D[i,j] * M₀[i,j]
    @tensor tmp_MV[i,j] = C[k] * F[i,j,k]
    @. Hartree = tmp_MV + newCᵨ/(Rmax-Rmin) * M₀
    discretization.cache.Cᵨ = newCᵨ
    nothing
end

#####################################################################
#                   Exchange Correlation Matrix
#####################################################################

function exchange_corr_matrix!(discretization::LDADiscretization, model, D)
    @unpack Vxc = discretization.cache
    ρ(x) = compute_density(discretization, D, x)
    weight(x) = vxc(model.exc, ρ(x))
    Vxc .= weight_mass_matrix(discretization.basis, weight)
    nothing
end

#####################################################################
#                             Energy
#####################################################################

function compute_energy!(discretization::LDADiscretization, solver::KhonShamSolver)
    compute_kinetic_energy!(discretization,solver)
    compute_coulomb_energy!(discretization,solver)
    compute_hartree_energy!(discretization,solver)
    if isthereExchangeCorrelation(solver.model)
        compute_exchangecorrelation_energy!(discretization,solver)
    end
    compute_total_energy!(discretization,solver)
end

function compute_total_energy!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack Rmax = discretization
    @unpack B, C, Cᵨ, Vxc = discretization.cache
    @unpack n, ϵ, D = solver
    @tensor energy = n[l,n] * ϵ[l,n] 
    if isthereExchangeCorrelation(solver.model)
        @tensor energy_correction = Vxc[i,j] * D[i,j]
        solver.energy = energy - discretization.elT(0.5) * (dot(B,C) + Cᵨ^2/Rmax) + solver.energy_exc - energy_correction
    else
        solver.energy = energy - discretization.elT(0.5) * (dot(B,C) + Cᵨ^2/Rmax)
    end
    nothing
end

function compute_kinetic_energy!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack A, M₋₂ = discretization.cache
    @unpack U, n = solver
    @unpack lₕ, Nₕ, elT  = discretization
    solver.energy_kin = zero(solver.energy_kin)
    @inbounds for l ∈ 1:lₕ+1   
        @inbounds for k ∈ 1:Nₕ
            if !iszero(n[l,k])
                solver.energy_kin += n[l,k] * elT(0.5) * U[l,:,k]' * (A + l*(l+1)*M₋₂) * U[l,:,k]
            end
        end
    end
    nothing
end

function compute_coulomb_energy!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack M₋₁ = discretization.cache
    @unpack U, n = solver
    @unpack lₕ, Nₕ  = discretization
    solver.energy_cou = zero(solver.energy_cou)
    @inbounds for l ∈ 1:lₕ+1   
        @inbounds for k ∈ 1:Nₕ
            if !iszero(n[l,k])
                solver.energy_cou -= solver.model.z * n[l,k] * U[l,:,k]' * M₋₁ * U[l,:,k]
            end
        end
    end
    nothing
end

function compute_hartree_energy!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack Rmax, elT = discretization
    @unpack B, C, Cᵨ = discretization.cache
    solver.energy_har = elT(0.5) * (dot(B,C) + Cᵨ^2/Rmax)
    nothing
end

function compute_exchangecorrelation_energy!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack D = solver
    @unpack Rmax = discretization
    ρ(x) = compute_density(discretization, D, x)
    f(x,p) = exc(solver.model.exc, ρ(x)) * x^2
    prob = IntegralProblem(f, (zero(Rmax),Rmax))
    solver.energy_exc = 4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
    nothing
end

#####################################################################
#                             Density
#####################################################################

function density_matrix!(discretization::LDADiscretization)
    @unpack tmp_Dstar, tmp_U, tmp_n = discretization.tmp_cache
    @unpack lₕ, Nₕ  = discretization
    @inbounds for l ∈ 1:lₕ+1   # potentiellement, il faut rajouter un lmin
        @inbounds for k ∈ 1 :Nₕ
            if !iszero(tmp_n[l,k])
                @inbounds for i ∈ 1:Nₕ
                    val = tmp_n[l,k] * tmp_U[l,i,k] 
                    @inbounds @simd for j ∈ 1:i
                        tmp_Dstar[i,j] += val * tmp_U[l,j,k]
                    end
                end
            end
        end
    end
    @inbounds for i in 1:Nₕ
        @inbounds @simd for j in 1:i-1
            tmp_Dstar[j,i] = tmp_Dstar[i,j]
        end
    end
    nothing
end

function compute_density(discretization::LDADiscretization, D, x)
    @unpack basis = discretization
    newT = promote_type(eltype(basis), typeof(x))
    val = zero(newT)
    eval_basis = zeros(newT, basis.size)
    @inbounds for i ∈ eachindex(basis)
        eval_basis[i] = basis(i,x)
    end
    val = (eval_basis)' * D * eval_basis
    return val* 1/4π * 1/(x^2)
end
