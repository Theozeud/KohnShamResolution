#####################################################################
#                          LDA Cache
#####################################################################

mutable struct LDAMatrices{T <: Real}
    # FEM MATRICES
    A::Matrix{T}                    # Matrix of Qᵢ'Qⱼ'
    M₀::Matrix{T}                   # Matrix of QᵢQⱼ
    M₋₁::Matrix{T}                  # Matrix of 1/x QᵢQⱼ
    M₋₂::Matrix{T}                  # Matrix of 1/x² QᵢQⱼ
    F::Array{T,3}                   # Tensor of 1/x QᵢQⱼQₖ
    # MATRICES COMPOSING THE HAMILTONIAN
    H::Array{T,3}                   # Hamiltonian
    Kin::Array{T,3}                 # Kinetic Matrix
    Coulomb::Matrix{T}              # Coulomb Matrix
    Hfix::Array{T,3}                # (Kinetic + Coulomb) Matrix        
    Hartree::Matrix{T}              # Hartree Matrix 
    Vxc::Matrix{T}                  # Echange-correlation Matrix
end

mutable struct LDACache{T <: Real}
    tmp_MV::Matrix{T}               # Store the contraction F:C
    tmp_B::Vector{T}                # Matrix of 4πρQᵢ  
    tmp_C::Vector{T}                # Matrix for V(ρ) - solution of the Gauss Electrostatic law
end


function create_cache_lda(lₕ::Int, Nₕ::Int, T)
    # FEM MATRICES
    A           = zeros(T, Nₕ, Nₕ) 
    M₀          = zeros(T, Nₕ, Nₕ)
    M₋₁         = zeros(T, Nₕ, Nₕ)
    M₋₂         = zeros(T, Nₕ, Nₕ)
    F           = zeros(T, Nₕ, Nₕ, Nₕ)

    # MATRICES COMPOSING THE HAMILTONIAN
    H           = zeros(T, lₕ+1, Nₕ, Nₕ)
    Kin         = zeros(T, lₕ+1, Nₕ, Nₕ)
    Coulomb     = zeros(T, Nₕ, Nₕ)
    Hfix        = zeros(T, lₕ+1, Nₕ, Nₕ)
    Hartree     = zeros(T, Nₕ, Nₕ)
    Vxc         = zeros(T, Nₕ, Nₕ)

    # Initialization of array for temporary stockage of computations 
    tmp_MV          = zeros(T, Nₕ, Nₕ)  
    tmp_B           = zeros(T, Nₕ)
    tmp_C           = zeros(T, Nₕ)

    LDAMatrices{T}(A, M₀, M₋₁, M₋₂, F, H, Kin, Coulomb, Hfix, Hartree, Vxc),  
    LDACache{T}(tmp_MV, tmp_B, tmp_C)
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
    matrices::LDAMatrices{T2}
    cache::LDACache{T2}
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

Base.eltype(discretization::LDADiscretization) = discretization.elT
dim(discretization::LDADiscretization) = discretization.Nₕ * (discretization.lₕ + 1)

#####################################################################
#                          Init Cache
#####################################################################

function init_cache!(discretization::LDADiscretization, model::AbstractDFTModel, hartree::Real)

    @unpack lₕ, basis, matrices  = discretization
    @unpack A, M₀, M₋₁, M₋₂, F, Kin, Coulomb, Hfix = matrices

    # CREATION OF FEM MATRICES
    fill_stiffness_matrix!(basis, A)
    fill_mass_matrix!(basis, M₀)
    fill_weight_mass_matrix!(basis, -1, M₋₁)
    lₕ == 0 || fill_weight_mass_matrix!(basis, -2, M₋₂)
    iszero(hartree) || fill_weight_mass_3tensor!(basis, Monomial(-1), F)

    # CREATION OF THE FIX PART OF THE HAMILTONIAN 
    kinetic_matrix!(discretization)
    coulomb_matrix!(discretization, model)
    for l ∈ 1:lₕ+1
        @views vHfix = Hfix[l,:,:]
        @views vKin = Kin[l,:,:]
        @. vHfix = vKin + Coulomb
    end

    nothing
end


#####################################################################
#                          Initialization
#####################################################################

init_density(kd::LDADiscretization)                 = zeros(kd.elT, kd.Nₕ, kd.Nₕ)  
init_orbitals(kd::LDADiscretization)                = zeros(kd.elT, kd.lₕ+1, kd.Nₕ, kd.Nₕ)
init_orbitals_energy(kd::LDADiscretization)         = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)
init_occupation_number(kd::LDADiscretization)       = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)
init_density_matrix(kd::LDADiscretization)          = BlockDiagonal([zeros(kd.elT, kd.Nₕ, kd.Nₕ) for i ∈ 1:kd.lₕ+1])

function init_energies(kd::LDADiscretization, model::KohnShamExtended)
    @unpack elT = kd
    d = Dict(   :Etot => zero(elT),                                     # Total energy 
                :Ekin => zero(elT),                                     # Kinetic energy
                :Ecou => zero(elT),                                     # Coulomb energy
                :Ehar => zero(elT))                                     # Hartree energy
    !isthereExchangeCorrelation(model) ||  d[:Eexc] = zero(elT)         # Exchange-correlation energy
    d
end


#####################################################################
#               Find Orbital : Solve the eigen problems
#####################################################################

function prepare_eigenvalue_problem!(   discretization::LDADiscretization, 
                                        model::KohnShamExtended, 
                                        D::AbstractMatrix{<:Real}, 
                                        hartree::Real = true)

    @unpack Hartree, Vxc = discretization.matrices

    # COMPUTE HARTREE MATRIX
    iszero(hartree) || hartree_matrix!(discretization, D, hartree)

    # COMPUTE EXCHANGE CORRELATION MATRIX
    !isthereExchangeCorrelation(model) || 
                    exchange_corr_matrix!(discretization, model, D)
    
    # BUILD THE HAMILTONIAN OF THE lᵗʰ SECTION
    @threads for l ∈ 0:lₕ
        @views vH = H[l+1,:,:]
        @views vHfix = Hfix[l+1,:,:]
        @. vH = vHfix + Vxc + Hartre
    end
    nothing
end

function find_orbital!( discretization::LDADiscretization, 
                        U::AbstractArray{<:Real}, 
                        ϵ::AbstractMatrix{<:Real})

    @unpack lₕ, matrices = discretization
    @unpack M₀, H = matrices

    # SOLVE THE GENERALIZED EIGENVALUE PROBLEM FOR EACH SECTION l
    @threads for l ∈ 0:lₕ
        @views vH = H[l+1,:,:]
        ϵ[l+1,:], U[l+1,:,:] = solve_generalized_eigenvalue_problem(vH, M₀)        
    end
end

#####################################################################
#               Normamisation of eigenvector
#####################################################################

function normalization!(discretization::LDADiscretization, 
                        U::AbstractArray{<:Real}, 
                        n::AbstractMatrix{<:Real})
    @unpack M₀ = discretization.matrices
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

function normalization!(discretization::LDADiscretization, U::AbstractArray{<:Real}, l::Int, k::Int)
    @unpack M₀ = discretization.matrices
    @views Ulk = U[l+1,:,k]
    coeff = sqrt(Ulk'*M₀*Ulk)
    Ulk .= Ulk .* 1.0/coeff
    nothing
end

#####################################################################
#                          Kinetic Matrix
#####################################################################

function kinetic_matrix!(discretization::LDADiscretization)
    @unpack A, M₋₂, Kin = discretization.matrices
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
    @unpack M₋₁, Coulomb = discretization.matrices
    Coulomb .= - model.z .* M₋₁
    nothing
end

#####################################################################
#                          Hartree Matrix
#####################################################################

function hartree_matrix!(discretization::LDADiscretization, D::AbstractMatrix{<:Real}, coeff::Real = true)
    @unpack Rmax, matrices, cache = discretization
    @unpack A, M₀, F, Hartree = matrices
    @unpack tmp_MV, tmp_B, tmp_C = cache
    @tensor tmp_B[m] = D[i,j] * F[i,j,m]
    tmp_C .= A\B
    @tensor newCrho = D[i,j] * M₀[i,j]
    @tensor tmp_MV[i,j] = tmp_C[k] * F[i,j,k]
    @. Hartree = tmp_MV + newCrho/Rmax * M₀
    @. Hartree .*= coeff
    nothing
end

#####################################################################
#                   Exchange Correlation Matrix
#####################################################################

function exchange_corr_matrix!( discretization::LDADiscretization, 
                                model::KohnShamExtended, 
                                D::AbstractMatrix{<:Real})
    @unpack matrices, basis = discretization
    @unpack Vxc = matrices
    ρ(x) = compute_density(discretization, D, x)
    weight(x) = vxc(model.exc, ρ(x))
    fill_weight_mass_matrix!(basis, weight, Vxc)
    nothing
end


#####################################################################
#                         TOTAL ENERGY
#####################################################################

function compute_total_energy!( discretization::LDADiscretization, 
                                model::KohnShamExtended,
                                D::AbstractMatrix{<:Real}, 
                                n::AbstractMatrix{<:Real},
                                ϵ::AbstractMatrix{<:Real})
    @unpack Rmax, matrices = discretization
    @unpack Vxc = matrices
    @tensor energy = n[l,n] * ϵ[l,n] 
    if isthereExchangeCorrelation(model)
        @tensor energy_correction = Vxc[i,j] * D[i,j]
        energy_exc = compute_exchangecorrelation_energy(discretization, model, D)
        energy_har = compute_hartree_energy(discretization, D)
        return energy - energy_har + energy_exc - energy_correction
    else
        energy_har = compute_hartree_energy(discretization, D)
        return energy = energy - energy_har
    end
    nothing
end

#####################################################################
#                        KINETIC ENERGY
#####################################################################

function compute_kinetic_energy(discretization::LDADiscretization, 
                                U::AbstractArray{<:Real}, 
                                n::AbstractArray{<:Real})
    @unpack lₕ, Nₕ, elT  = discretization
    @unpack Kin = discretization.matrices
    energy_kin = zero(elT)
    @inbounds for l ∈ 1:lₕ+1 
        @views vKin = Kin[l,:,:]  
        @inbounds for k ∈ 1:Nₕ
            if !iszero(n[l,k])
                @views Ulk = U[l,:,k]
                energy_kin += n[l,k] * Ulk' * vKin * Ulk
            end
        end
    end
    return energy_kin
end

#####################################################################
#                        COULOMB ENERGY
#####################################################################

function compute_coulomb_energy(discretization::LDADiscretization, U::AbstractArray{<:Real}, n::AbstractArray{<:Real})
    @unpack lₕ, Nₕ, elT  = discretization
    @unpack Coulomb = discretization.matrices
    energy_cou = zero(elT)
    @inbounds for l ∈ 1:lₕ+1   
        @inbounds for k ∈ 1:Nₕ
            if !iszero(n[l,k])
                @views Ulk = U[l,:,k]
                energy_cou +=  n[l,k] * Ulk' * Coulomb * Ulk
            end
        end
    end
    return energy_cou
end

#####################################################################
#                        HARTREE ENERGY
#####################################################################

function compute_hartree_energy(discretization::LDADiscretization, D::AbstractMatrix{<:Real})
    @unpack Rmax, elT, matrices, cache = discretization
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    @tensor tmp_B[m] = D[i,j] * F[i,j,m]
    tmp_C .= A\tmp_B
    @tensor Crho = D[i,j] * M₀[i,j]
    return elT(0.5) * (dot(tmp_B,tmp_C) + Crho^2/Rmax)
end

function compute_hartree_mix_energy(discretization::LDADiscretization, 
                                    D0::AbstractMatrix{<:Real}, 
                                    D1::AbstractMatrix{<:Real})
    @unpack Rmax, elT, matrices, cache = discretization
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    @tensor tmp_B[m] = D0[i,j] * F[i,j,m]
    tmp_C .= A\tmp_B
    @tensor tmp_B[m] = D1[i,j] * F[i,j,m]
    @tensor Crho0 = D0[i,j] * M₀[i,j]
    @tensor Crho1 = D1[i,j] * M₀[i,j]
    return elT(0.5) * (dot(tmp_B,tmp_C) + Crho0*Crho1/Rmax)
end

#####################################################################
#                  EXCHANGE CORRELATION ENERGY
#####################################################################

function compute_exchangecorrelation_energy(discretization::LDADiscretization, 
                                            model::KohnShamExtended, 
                                            D::AbstractMatrix{<:Real})
    @unpack Rmax = discretization
    ρ(x) = compute_density(discretization, D, x)
    f(x,p) = exc(model.exc, ρ(x)) * x^2
    prob = IntegralProblem(f, (zero(Rmax),Rmax))
    4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
end

#####################################################################
#                             Density
#####################################################################

function density!(  discretization::LDADiscretization, 
                    U::AbstractArray{<:Real}, 
                    n::AbstractMatrix{<:Real}, 
                    D::AbstractMatrix{<:Real})
    @unpack lₕ, Nₕ, elT  = discretization
    fill!(D, zero(elT))
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

function density!(  discretization::LDADiscretization, 
                    Γ::BlockDiagonal{<:Real, <:AbstractMatrix{<:Real}},
                    D::AbstractMatrix{<:Real})
    @unpack lₕ, Nₕ, elT  = discretization
    fill!(D, zero(elT))
    @inbounds @simd for l ∈ 1:lₕ+1
        @views Γl = blocks(Γ)[l]
        @inbounds for i ∈ 1:Nₕ
            @inbounds for j ∈ 1:i
                D[i,j] += (2*l + 1)*Γl[i,j]
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

function compute_density(discretization::LDADiscretization, D::AbstractMatrix{<:Real}, x::Real)
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


#####################################################################
#                          Density Matrix
#####################################################################

function density_matrix!(   discretization::LDADiscretization, 
                            U::AbstractArray{<:Real}, 
                            n::AbstractMatrix{<:Real}, 
                            Γ::BlockDiagonal{<:Real, <:AbstractMatrix{<:Real}})
    @unpack lₕ, Nₕ, elT  = discretization
    @inbounds for l ∈ 1:lₕ+1 
        @views Γl = blocks(Γ)[l]
        fill!(Γl, zero(elT))
        @inbounds for k ∈ 1:Nₕ
            if !iszero(n[l,k])
                @inbounds for i ∈ 1:Nₕ
                    val = n[l,k]/(2*l-1) * U[l,i,k] 
                    @inbounds @simd for j ∈ 1:i
                        Γl[i,j] += val * U[l,j,k]
                    end
                end
            end
        end
        @inbounds for i in 1:Nₕ
            @inbounds @simd for j in 1:i-1
                Γl[j,i] = Γl[i,j]
            end
        end
    end
    nothing
end