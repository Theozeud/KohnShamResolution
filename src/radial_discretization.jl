#####################################################################
#                          Radial Cache
#####################################################################

mutable struct RadialCache
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
    Exc             # Matrix VF of Exc
    Energy          # Total Energy
end

mutable struct Radial_tmp_Cache
    tmp_H           # Store Hₗ
    tmp_D           # Store ρ
    tmp_Dstar       # Store ρ*
    tmp_U           # Store the eigenvector
    tmp_ϵ           # Store the eigenvalue
    tmp_n           # Store the occupation number
    tmp_MV          # Store the contraction F:C
    tmp_index_sort  #
    
end

function create_cache(lₕ, Nₕ, T, lmin)
    A       = zeros(T, Nₕ, Nₕ) 
    M₀      = zeros(T, Nₕ, Nₕ)
    M₋₁     = zeros(T, Nₕ, Nₕ)
    M₋₂     = zeros(T, Nₕ, Nₕ)
    F       = zeros(T, Nₕ, Nₕ, Nₕ)
    B       = zeros(T, Nₕ)
    C       = zeros(T, Nₕ)
    Cᵨ      = zero(T)
    Kin     = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    Coulomb = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    Hfix    = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    Hartree = zeros(T, Nₕ, Nₕ)
    Exc     = zeros(T, Nₕ, Nₕ)
    Energy  = zero(T)

    # Initialization of array for temporary stockage of computations
    tmp_H           = zeros(T, Nₕ, Nₕ)
    tmp_D           = zeros(T, Nₕ, Nₕ) 
    tmp_Dstar       = zeros(T, Nₕ, Nₕ) 
    tmp_U           = zeros(T, lₕ+1 - lmin, Nₕ, Nₕ)
    tmp_MV          = zeros(T, Nₕ, Nₕ)
    tmp_ϵ           = zeros(T, lₕ+1 - lmin, Nₕ)
    tmp_index_sort  = zeros(Int, Nₕ*(lₕ+1 - lmin))
    tmp_n           = zeros(T, lₕ+1 - lmin, Nₕ)   
 
    RadialCache(A, M₀, M₋₁, M₋₂, F, B, C, Cᵨ, Kin, Coulomb, Hfix, Hartree, Exc, Energy),  Radial_tmp_Cache(tmp_H, tmp_D, tmp_Dstar, tmp_U,  tmp_ϵ, tmp_n, tmp_MV, tmp_index_sort)
end


#####################################################################
#                          Radial Discretization
#####################################################################
abstract type KohnShamDiscretization end

struct KohnShamRadialDiscretization{T} <: KohnShamDiscretization
    lmin::Int
    lₕ::Int
    Nₕ::Int
    basis::Basis
    mesh::Mesh{T}
    Rmin::T
    Rmax::T
    elT::Type
    cache::RadialCache
    tmp_cache::Radial_tmp_Cache
    function KohnShamRadialDiscretization(lₕ::Int, basis::Basis, mesh::Mesh; lmin = 0)
        elT = bottom_type(basis)
        Nₕ = length(basis)
        @assert lmin ≤ lₕ
        new{eltype(mesh)}(lmin, lₕ, Nₕ, basis, mesh, first(mesh), last(mesh), elT, create_cache(lₕ, Nₕ, elT, lmin)...)
    end
end

#####################################################################
#                          Init Cache
#####################################################################

function init_cache!(discretization::KohnShamRadialDiscretization, model::AbstractDFTModel, hartree)

    @unpack lₕ, basis  = discretization
    @unpack A, M₀, M₋₁, M₋₂, F, Kin, Coulomb, Hfix = discretization.cache

    # Creation of the base matrices
    deriv_basis = deriv(basis)
    A   .= mass_matrix(deriv_basis)
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

init_density_matrix(kd::KohnShamRadialDiscretization)        = zeros(kd.elT, kd.Nₕ, kd.Nₕ)  
init_coeffs_discretization(kd::KohnShamRadialDiscretization) = zeros(kd.elT, kd.lₕ+1, kd.Nₕ, kd.Nₕ)
init_energy(kd::KohnShamRadialDiscretization)                = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)
init_occupation(kd::KohnShamRadialDiscretization)            = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)

#####################################################################
#               Find Orbital : Solve the eigen problems
#####################################################################

function find_orbital!(discretization::KohnShamRadialDiscretization, solver::KhonShamSolver)

    @unpack lmin, lₕ = discretization
    @unpack M₀, Hfix, Hartree, Exc = discretization.cache
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
    for l ∈ lmin:lₕ
        # building the hamiltonian of the lᵗʰ section
        @. tmp_H = Hfix[l+1-lmin,:,:] + Exc + Hartree
        # solving
        tmp_ϵ[l+1-lmin,:], tmp_U[l+1-lmin,:,:] = solve_generalized_eigenvalue_problem(tmp_H, M₀)
    end
end

#####################################################################
#                          Kinetic Matrix
#####################################################################

function kinetic_matrix!(discretization::KohnShamRadialDiscretization)
    @unpack A, M₋₂, Kin = discretization.cache
    for l ∈ discretization.lmin:discretization.lₕ
        @. Kin[l+1-discretization.lmin,:,:] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end

#####################################################################
#                          Coulomb Matrix
#####################################################################

function coulomb_matrix!(discretization::KohnShamRadialDiscretization, model)
    @unpack M₋₁, Coulomb = discretization.cache
    for l ∈ discretization.lmin:discretization.lₕ
        Coulomb[l+1-discretization.lmin,:,:] .= - model.z .* M₋₁
    end 
    nothing
end

#####################################################################
#                          Hartree Matrix
#####################################################################

function hartree_matrix!(discretization::KohnShamRadialDiscretization, D)
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

function exchange_corr_matrix!(discretization::KohnShamRadialDiscretization, model, D)
    @unpack Exc = discretization.cache
    ρ(x) = eval_density_as_function(discretization, D, x)
    Exc .= weight_mass_matrix(discretization.basis, model.exc.vxc∘ρ)
    nothing
end

#####################################################################
#                             Energy
#####################################################################

function energy(discretization::KohnShamRadialDiscretization)
    @unpack Rmin, Rmax = discretization
    @unpack Energy, B, C, Cᵨ = discretization.cache
    @unpack tmp_n, tmp_ϵ = discretization.tmp_cache
    @tensor Energy = tmp_n[l,n] * tmp_ϵ[l,n] 
    discretization.cache.Energy = Energy - discretization.elT(0.5) * (dot(B,C) + Cᵨ^2/(Rmax-Rmin))
    nothing
end

#####################################################################
#                         Eigenvector
#####################################################################

function build_eigenvector(kd::KohnShamRadialDiscretization, U; Index = CartesianIndices((1:lₕ+1, 1:Nₕ)))
    @unpack lₕ, Nₕ, basis, mesh = kd
    # First computation is done separately to well instantiate the array of eigenvector
    Ifirst = first(Index)
    lfirst = Ifirst[1]
    kfirst = Ifirst[2]
    eiglkfirst = build_on_basis(basis, U[lfirst,:,kfirst])
    eigenvectors = [1/sqrt(4π)* eiglkfirst / normL2(eiglkfirst, mesh) * Monomial(-1)]
    for I ∈ Index[begin+1:end]
        l = I[1]
        k = I[2]
        eiglk = build_on_basis(basis, U[l,:,k]) 
        push!(eigenvectors,  1/sqrt(4π)* eiglk / normL2(eiglk, mesh) * Monomial(-1)) 
    end
    eigenvectors
end

#####################################################################
#                             Density
#####################################################################

function density_matrix!(discretization::KohnShamRadialDiscretization)
    @unpack M₀ = discretization.cache
    @unpack tmp_Dstar, tmp_U, tmp_n = discretization.tmp_cache
    @unpack lₕ, Nₕ  = discretization
    @inbounds for l ∈ 1:lₕ+1 
        @inbounds for k ∈ 1 :Nₕ
            if !iszero(tmp_n[l,k])
                normalization = sum([tmp_U[l,i,k] * tmp_U[l,j,k] * M₀[i,j] for i∈1:Nₕ for j∈1:Nₕ])
                @inbounds for i ∈ 1:Nₕ
                    val = tmp_n[l,k] * tmp_U[l,i,k] * 1/normalization
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

function build_density!(discretization::KohnShamRadialDiscretization, D)
    @unpack basis = discretization
    ρ = Monomial(0,0)
    Basis = [build_basis(basis, i) for i∈axes(D,1)]
    for i ∈ axes(D,1)
        for j ∈ axes(D,2)
            ρ += D[i,j]  * Basis[i] * Basis[j]
        end
    end
    ρ* 1/4π * Monomial(-2)
end

function eval_density_as_function(discretization::KohnShamRadialDiscretization, D, x)
    @unpack basis, mesh = discretization
    if iszero(x)
        x+= mesh[2]/10
    end
    val = 0
    for i ∈ axes(D,1)
        vx = eval_basis(basis, i, x)
        for j ∈ axes(D,2)
            vy = eval_basis(basis, i, x)
            val += D[i,j]  * vx * vy
        end
    end
    return val* 1/4π * 1/x^2
end



