#####################################################################
#                          Radial Cache
#####################################################################

struct RadialCache
    A
    M₀
    M₋₁
    M₋₂
    F
    B
    Kin
    Coulomb
    Hfix
    C
    Cᵨ
    Hartree
    Exc
    Energy
    tmp_H
    tmp_D
    tmp_Dstar
    tmp_U
    tmp_MV
    tmp_ϵ
    tmp_index_sort
    tmp_n    
end

function create_cache(lₕ, Nₕ, T)

    A   = zeros(T, Nₕ, Nₕ) 
    M₀  = zeros(T, Nₕ, Nₕ)
    M₋₁ =  zeros(T, Nₕ, Nₕ)
    M₋₂ =  zeros(T, Nₕ, Nₕ)
    F   =  zeros(T, Nₕ, Nₕ, Nₕ)
    B   =  zeros(T, Nₕ)
    Kin =  zeros(T, lₕ+1, Nₕ, Nₕ)
    Coulomb = zeros(T, lₕ+1, Nₕ, Nₕ)
    Hfix    = zeros(T, lₕ+1, Nₕ, Nₕ)
    C     = zeros(T, Nₕ)
    Cᵨ    = zero(T)
    Hartree = zeros(T, Nₕ, Nₕ)
    Exc     = zeros(T, Nₕ, Nₕ)
    Energy  = zero(T)

    # Initialization of array for temporary stockage of computations
    tmp_H           = zeros(T, Nₕ, Nₕ)
    tmp_D           = zeros(T, Nₕ, Nₕ) #zero_piecewiselaurantpolynomial(mesh, T)
    tmp_Dstar       = zeros(T, Nₕ, Nₕ) #zero_piecewiselaurantpolynomial(mesh, T)
    tmp_U           = zeros(T, lₕ+1, Nₕ, Nₕ)
    tmp_MV          = zeros(T, Nₕ, Nₕ)
    tmp_ϵ           = zeros(T, lₕ+1, Nₕ)
    tmp_index_sort  = zeros(Int, Nₕ*(lₕ+1))
    tmp_n           = zeros(T, lₕ+1, Nₕ)   
 
    RadialCache(A, M₀, M₋₁, M₋₂, F, B, Kin, Coulomb, Hfix, C, Cᵨ, Hartree, Exc, Energy, tmp_H, tmp_D, tmp_Dstar, tmp_U, tmp_MV, tmp_ϵ, tmp_index_sort, tmp_n)
end


#####################################################################
#                          Radial Discretization
#####################################################################


struct KohnShamRadialDiscretization{T} <: KohnShamDiscretization
    lₕ::Int
    Nₕ::Int
    basis::Basis
    mesh::OneDMesh{T}
    Rmin::T
    Rmax::T
    elT::Type
    cache::RadialCache
    function KohnShamRadialDiscretization(lₕ::Int, basis::Basis, mesh::OneDMesh)
        elT = bottom_type(basis)
        Nₕ = length(basis)
        new{eltype(mesh)}(lₕ, Nₕ, basis, mesh, first(mesh), last(mesh), elT, create_cache(lₕ, Nₕ, elT))
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
    if hartree
        F .= weight_mass_3tensor(basis, Monomial(-1))
    end

    nothing
end

#init_density_matrix(kd::KohnShamRadialDiscretization)        = zero_piecewiselaurantpolynomial(kd.mesh, kd.elT), zero_piecewiselaurantpolynomial(kd.mesh, kd.elT)
init_density_matrix(kd::KohnShamRadialDiscretization)        = zeros(kd.elT, kd.Nₕ, kd.Nₕ)  
init_coeffs_discretization(kd::KohnShamRadialDiscretization) = zeros(kd.elT, kd.lₕ+1, kd.Nₕ, kd.Nₕ)
init_energy(kd::KohnShamRadialDiscretization)                = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)
init_occupation(kd::KohnShamRadialDiscretization)            = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)

#####################################################################
#                          Kinetic Energy
#####################################################################

function kinetic_matrix!(discretization::KohnShamRadialDiscretization)
    @unpack A, M₋₂, Kin = discretization.cache
    for l ∈ 0:discretization.lₕ
        @. Kin[l+1,:,:] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end

#####################################################################
#                          Coulomb Energy
#####################################################################

function coulomb_matrix!(discretization::KohnShamRadialDiscretization, model)
    @unpack M₋₁, Coulomb = discretization.cache
    for l ∈ 0:discretization.lₕ
        Coulomb[l+1,:,:] .= - charge(model) .* M₋₁
    end 
    nothing
end

#####################################################################
#                          Hartree Energy
#####################################################################

function hartree_matrix!(discretization::KohnShamRadialDiscretization, ρ, opt)
    if opt == :integral
        #build_hartree_integral!(kd::KohnShamRadialDiscretization, Hartree, ρ)
    elseif opt == :pde
        hartree_matrix_pde!(discretization, ρ)
    end
    discretization.cache.tmp_Hartree
end

function hartree_matrix_integral!(kd::KohnShamRadialDiscretization, ρ)
    @unpack basis, Rmin, Rmax = kd
    potential = nothing
    fill_weight_mass_matrix!(basis, potential, Hartree)
    Hartree
end

function hartree_matrix_pde_old!(discretization::KohnShamRadialDiscretization, ρ)
    @unpack A, M₀, Hartree = discretization.cache
    @unpack basis, Rmin, Rmax = discretization
    rho(x) = ρ(x)
    f(x) = 4π * rho(x) * x
    F =  weight_mass_vector(basis, f)
    coeff = A\F
    Cᵨ = 4π * integrate(ρ * Monomial(2), Rmin, Rmax)
    MV  = vectorweight_mass_matrix(basis, coeff, Monomial(-1)) 
    @. Hartree = MV + Cᵨ/(Rmax-Rmin) * M₀
    nothing
end

function hartree_matrix_pde!(discretization::KohnShamRadialDiscretization, D)
    @unpack A, M₀, F, B, C, Cᵨ, Hartree, tmp_MV = discretization.cache
    @unpack basis, Rmin, Rmax = discretization
    @tensor B[m] = D[i,j] * F[i,j,m]
    C .= A\B
    @tensor Cᵨ = D[i,j] * M₀[i,j]
    @tensor tmp_MV[i,j] = C[k] * F[i,j,k]
    @. Hartree = tmp_MV + Cᵨ/(Rmax-Rmin) * M₀
    nothing
end

#####################################################################
#                         Exchange Correlation
#####################################################################

function exchange_corr_matrix!(discretization::KohnShamRadialDiscretization, model)
    @unpack Exc = discretization.cache
    Exc .= weight_mass_matrix(basis, model.exc)
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
    @unpack M₀, tmp_Dstar, tmp_U, tmp_n = discretization.cache
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

function build_density!(discretization::KohnShamRadialDiscretization)
    @unpack tmp_Dstar, tmp_U, tmp_n = discretization.cache
    @unpack lₕ, Nₕ, basis, mesh = discretization
    for l ∈ 1:lₕ+1
        for k ∈ 1:Nₕ
            if tmp_n[l,k] != 0
                eigen_vector = build_on_basis(basis, tmp_U[l,:,k]) 
                eigen_vector = eigen_vector / normL2(eigen_vector, mesh) * Monomial(-1)
                tmp_Dstar += tmp_n[l,k] * eigen_vector * eigen_vector
            end
        end
    end
    tmp_Dstar *= 1/(4π)
end

function build_density2!(discretization::KohnShamRadialDiscretization)
    @unpack tmp_Dstar = discretization.cache
    ρ = Monomial(0,0)
    for i ∈ axes(D,1)
        for j ∈ axes(D,2)
            ρ += D[i,j]  * build_basis(basis, i) * build_basis(basis, j)
        end
    end
    ρ* 1/4π * Monomial(-2)
end


#####################################################################
#                             Energy
#####################################################################

function energy(discretization::KohnShamRadialDiscretization)
    @unpack Rmax = discretization
    @unpack Energy, tmp_n, tmp_ϵ, B, C, Cᵨ = discretization.cache
    @tensor Energy = tmp_n[l,n] * tmp_ϵ[l,n] 
    Energy - discretization.elT(0.5) * (dot(B,C) + Cᵨ^2/Rmax)
end

