#####################################################################
#                          LDA Cache
#####################################################################

mutable struct LDACache{T <: Real}
    A::Matrix{T}                    # Matrix of Qᵢ'Qⱼ'
    M₀::Matrix{T}                   # Matrix of QᵢQⱼ
    M₋₁::Matrix{T}                  # Matrix of 1/x QᵢQⱼ
    M₋₂::Matrix{T}                  # Matrix of 1/x^2 QᵢQⱼ
    F::Array{T}                     # Tensor of 1/x QᵢQⱼQₖ
    B::Vector{T}                    # Matrix of 4πρQᵢ
    C::Vector{T}                    # Matrix for V(ρ) - solution of the Gauss Electrostatic law
    Cprev::Vector{T}                # Matrix for V(ρ) - solution of the Gauss Electrostatic law at previous time
    Cᵨ::T                           # Number equal to total charge
    Cᵨprev::T                       # Number equal to total charge at previous time
    Kin::Array{T}                  # Matrix VF of Kinetic 
    Coulomb::Matrix{T}              # Matrix VF of Coulomb
    Hfix::Array{T}                 # Part of Hamilotnian not needing to be recomputed (Kinetic + Colombial)         
    Hartree::Matrix{T}              # Matrix VF of hartree 
    Vxc::Matrix{T}                  # Matrix VF of Echange-correlation
end

mutable struct LDA_tmp_Cache{T <: Real}
    tmp_H::Array{T}                 # Store Hₗ
    tmp_MV::Matrix{T}               # Store the contraction F:C  
    tmp_index_sort::Vector{Int}
end

function create_cache_lda(lₕ::Int, Nₕ::Int, T)
    A           = zeros(T, Nₕ, Nₕ) 
    M₀          = zeros(T, Nₕ, Nₕ)
    M₋₁         = zeros(T, Nₕ, Nₕ)
    M₋₂         = zeros(T, Nₕ, Nₕ)
    F           = zeros(T, Nₕ, Nₕ, Nₕ)
    B           = zeros(T, Nₕ)
    C           = zeros(T, Nₕ)
    Cprev       = zeros(T, Nₕ)
    Cᵨ          = zero(T)
    Cᵨprev      = zero(T)
    Kin         = zeros(T, lₕ+1, Nₕ, Nₕ)
    Coulomb     = zeros(T, Nₕ, Nₕ)
    Hfix        = zeros(T, lₕ+1, Nₕ, Nₕ)
    Hartree     = zeros(T, Nₕ, Nₕ)
    Vxc         = zeros(T, Nₕ, Nₕ)

    # Initialization of array for temporary stockage of computations
    tmp_H           = zeros(T, lₕ+1, Nₕ, Nₕ) 
    tmp_MV          = zeros(T, Nₕ, Nₕ)  
    tmp_index_sort  = zeros(Int, Nₕ*(lₕ+1))
    LDACache{T}(A, M₀, M₋₁, M₋₂, F, B, C, Cprev, Cᵨ, Cᵨprev, Kin, Coulomb, Hfix, Hartree, Vxc),  
    LDA_tmp_Cache{T}(tmp_H, tmp_MV, tmp_index_sort)
end


#####################################################################
#                          LDA Discretization
#####################################################################


struct LDADiscretization{T <: Real, T2 <: Real, typeBasis <: Basis} <: KohnShamDiscretization
    lₕ::Int
    Nₕ::Int
    basis::typeBasis
    mesh::Mesh{T}
    Rmin::T
    Rmax::T
    elT::Type
    cache::LDACache{T2}
    tmp_cache::LDA_tmp_Cache{T2}
    function LDADiscretization(lₕ::Int, basis::Basis, mesh::Mesh)
        elT = try
             bottom_type(basis)
        catch
            eltype(basis)
        end
        Nₕ = length(basis)
        new{eltype(mesh), elT, typeof(basis)}(lₕ, Nₕ, basis, mesh, first(mesh), last(mesh), elT, create_cache_lda(lₕ, Nₕ, elT)...)
    end
end

#####################################################################
#                          Init Cache
#####################################################################

function init_cache!(discretization::LDADiscretization, model::AbstractDFTModel, hartree::Real)

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
    for l ∈ 1:lₕ+1
        @views vHfix = Hfix[l,:,:]
        @views vKin = Kin[l,:,:]
        @. vHfix = vKin + Coulomb
    end

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

    @unpack lₕ = discretization
    @unpack M₀, Hfix, Hartree, Vxc = discretization.cache
    @unpack tmp_H = discretization.tmp_cache
    @unpack Dprev, U, ϵ = solver
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
    @threads for l ∈ 0:lₕ
        # building the hamiltonian of the lᵗʰ section
        @views vtmpH = tmp_H[l+1,:,:]
        @views vHfix = Hfix[l+1,:,:]
        vtmpH .= vHfix + Vxc + Hartree
        # solving
        ϵ[l+1,:], U[l+1,:,:] = solve_generalized_eigenvalue_problem(vtmpH, M₀)        
    end
end

#####################################################################
#               Normamisation of eigenvector
#####################################################################

function normalization!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack M₀ = discretization.cache
    @unpack U, n = solver
    @unpack lₕ, Nₕ  = discretization
    @inbounds for k ∈ 1:Nₕ
        @inbounds for l ∈ 1:lₕ+1   
            if !iszero(n[l,k])
                @views Ulk = U[l,:,k]
                coeff = sqrt(Ulk'*M₀*Ulk)
                Ulk .= Ulk .* 1.0/coeff
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
    for l ∈ 0:discretization.lₕ
        @views vkin = Kin[l+1,:,:]
        @. vkin =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end

#####################################################################
#                          Coulomb Matrix
#####################################################################

function coulomb_matrix!(discretization::LDADiscretization, model::KohnShamExtended)
    @unpack M₋₁, Coulomb = discretization.cache
    Coulomb .= - model.z .* M₋₁
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

function exchange_corr_matrix!(discretization::LDADiscretization, model::KohnShamExtended, D)
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
    @unpack Kin = discretization.cache
    @unpack U, n = solver
    @unpack lₕ, Nₕ, elT  = discretization
    solver.energy_kin = zero(solver.energy_kin)
    @inbounds for l ∈ 1:lₕ+1 
        @views vKin = Kin[l,:,:]  
        @inbounds for k ∈ 1:Nₕ
            if !iszero(n[l,k])
                @views Ulk = U[l,:,k]
                solver.energy_kin += n[l,k] * Ulk' * vKin * Ulk
            end
        end
    end
    nothing
end

function compute_coulomb_energy!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack Coulomb = discretization.cache
    @unpack U, n = solver
    @unpack lₕ, Nₕ  = discretization
    solver.energy_cou = zero(solver.energy_cou)
    @inbounds for l ∈ 1:lₕ+1   
        @inbounds for k ∈ 1:Nₕ
            if !iszero(n[l,k])
                @views Ulk = U[l,:,k]
                solver.energy_cou -=  n[l,k] * Ulk' * Coulomb * Ulk
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

function compute_hartree_mix_energy(discretization::LDADiscretization, ::KhonShamSolver)
    @unpack Rmax, elT = discretization
    @unpack B, Cᵨ, Cᵨprev, Cprev = discretization.cache
    return elT(0.5) * (dot(B,Cprev) + Cᵨ*Cᵨprev/Rmax)
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

function density_matrix!(discretization::LDADiscretization, solver::KhonShamSolver)
    @unpack U, n, D = solver
    @unpack lₕ, Nₕ  = discretization
    @inbounds for k ∈ 1 :Nₕ
        @inbounds for l ∈ 1:lₕ+1   
            if !iszero(n[l,k])
                @inbounds for i ∈ 1:Nₕ
                    val = n[l,k] * U[l,i,k] 
                    @inbounds @simd for j ∈ 1:i
                        D[i,j] += val * U[l,j,k]
                    end
                end
            end
        end
    end
    @inbounds for i in 1:Nₕ
        @inbounds @simd for j in 1:i-1
            D[j,i] = D[i,j]
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
