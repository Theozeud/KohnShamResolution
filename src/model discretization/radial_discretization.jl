#####################################################################
#                          Radial Cache
#####################################################################

struct RadialCache
    A
    M₀
    M₋₁
    M₋₂
    Kin
    Coulomb
    Hfix
    tmp_H
    tmp_D
    tmp_Dstar
    tmp_U
    tmp_Hartree
    tmp_exc 
    tmp_ϵ
    tmp_index_sort
    tmp_n    
end

function create_cache(lₕ, Nₕ, mesh, T)

    A   = zeros(T, Nₕ, Nₕ) 
    M₀  = zeros(T, Nₕ, Nₕ)
    M₋₁ =  zeros(T, Nₕ, Nₕ)
    M₋₂ =  zeros(T, Nₕ, Nₕ)
    Kin =  zeros(T, lₕ+1, Nₕ, Nₕ)
    Coulomb =  zeros(T, lₕ+1, Nₕ, Nₕ)
    Hfix = zeros(T, lₕ+1, Nₕ, Nₕ) 

    # Initialization of array for temporary stockage of computations
    tmp_H           = zeros(T, Nₕ, Nₕ)
    tmp_D           = zero_piecewiselaurantpolynomial(mesh, T)
    tmp_Dstar       = zero_piecewiselaurantpolynomial(mesh, T)
    tmp_U           = zeros(T, lₕ+1, Nₕ, Nₕ)
    tmp_Hartree     = zeros(T, Nₕ, Nₕ)
    tmp_exc         = zeros(T, Nₕ, Nₕ)
    tmp_ϵ           = zeros(T, lₕ+1, Nₕ)
    tmp_index_sort  = zeros(Int, Nₕ*(lₕ+1))
    tmp_n           = zeros(T, lₕ+1, Nₕ)   
 
    RadialCache(A, M₀, M₋₁, M₋₂, Kin, Coulomb, Hfix, tmp_H, tmp_D, tmp_Dstar, tmp_U, tmp_Hartree, tmp_exc, tmp_ϵ, tmp_index_sort, tmp_n)
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
        new{eltype(mesh)}(lₕ, Nₕ, basis, mesh, first(mesh), last(mesh), elT, create_cache(lₕ, Nₕ, mesh, elT))
    end
end

#####################################################################
#                          Init Cache
#####################################################################

function init_cache!(discretization::KohnShamRadialDiscretization, model::AbstractDFTModel)

    @unpack lₕ, basis  = discretization
    @unpack A, M₀, M₋₁, M₋₂, Kin, Coulomb, Hfix = discretization.cache

    # Creation of the base matrices
    deriv_basis = deriv(basis)
    A   .= mass_matrix(deriv_basis)
    M₀  .= mass_matrix(basis)
    M₋₁ .= weight_mass_matrix(basis, -1)

    # Creation of the fix part of the hamiltonian   
    if lₕ ≠ 0
        M₋₂ .= weight_mass_matrix(basis, -2)
    end
    build_kinetic!(discretization)
    build_coulomb!(discretization, model)
    @. Hfix = Kin + Coulomb
    nothing
end

init_density_matrix(kd::KohnShamRadialDiscretization)        = zero_piecewiselaurantpolynomial(kd.mesh, kd.elT), zero_piecewiselaurantpolynomial(kd.mesh, kd.elT) 
init_coeffs_discretization(kd::KohnShamRadialDiscretization) = zeros(kd.elT, kd.lₕ+1, kd.Nₕ, kd.Nₕ)
init_energy(kd::KohnShamRadialDiscretization)                = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)
init_occupation(kd::KohnShamRadialDiscretization)            = zeros(kd.elT, kd.lₕ+1, kd.Nₕ)

#####################################################################
#                          Kinetic Energy
#####################################################################

function build_kinetic!(discretization::KohnShamRadialDiscretization)
    @unpack A, M₋₂, Kin = discretization.cache
    for l ∈ 0:discretization.lₕ
        @. Kin[l+1,:,:] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end


#####################################################################
#                          Coulomb Energy
#####################################################################

function build_coulomb!(discretization::KohnShamRadialDiscretization, model)
    @unpack M₋₁, Coulomb = discretization.cache
    for l ∈ 0:discretization.lₕ
        Coulomb[l+1,:,:] .= - charge(model) .* M₋₁
    end 
    nothing
end

#####################################################################
#                          Hartree Energy
#####################################################################

function build_hartree!(discretization::KohnShamRadialDiscretization, ρ, opt)
    if opt == :integral
        #build_hartree_integral!(kd::KohnShamRadialDiscretization, Hartree, ρ)
    elseif opt == :pde
        build_hartree_pde!(discretization, ρ)
    end
    discretization.cache.tmp_Hartree
end

function build_hartree_integral!(kd::KohnShamRadialDiscretization, ρ)
    @unpack basis, Rmin, Rmax = kd
    potential = nothing
    fill_weight_mass_matrix!(basis, potential, Hartree)
    Hartree
end

function build_hartree_pde!(discretization::KohnShamRadialDiscretization, ρ)
    @unpack A, M₀, tmp_Hartree = discretization.cache
    @unpack basis, Rmin, Rmax = kd
    rho(x) = ρ(x)
    f(x) = 4π * rho(x) * x
    F =  weight_mass_vector(basis, f)
    coeff = A\F
    Cᵨ = 4π * integrate(ρ * Monomial(2), Rmin, Rmax)
    MV  = vectorweight_mass_matrix(basis, coeff, Monomial(-1)) 
    @. tmp_Hartree = MV + Cᵨ/(Rmax-Rmin) * M₀
    nothing
end

#####################################################################
#                         Exchange Correlation
#####################################################################

function build_exchange_corr!(kd::KohnShamRadialDiscretization, exc_mat, ρ, exc::AbstractExchangeCorrelation; quad_method, quad_reltol, quad_abstol)
    @unpack Nₕ, basis, Rmin, Rmax = kd
    if !isthereExchangeCorrelation(exc)
        return exc_mat .= zeros(Nₕ,Nₕ)
    end
    for i ∈ eachindex(basis)
        for j ∈ eachindex(basis)
            if j<i
                exc_mat[i,j] = exc_mat[j,i]
            else
                exc_mat[i,j] = approximate_integral(x -> exc.vxc(ρ(x)) * basis[i](x) * basis[j](x), (Rmin, Rmax) ; method = quad_method, reltol = quad_reltol, abstol = quad_abstol)
            end
        end
    end
    exc_mat
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
    nothing
end

#####################################################################
#                             Energy
#####################################################################

function build_energy()


end

