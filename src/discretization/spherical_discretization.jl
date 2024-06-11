struct KohnShamSphericalDiscretization{T} <: KohnShamDiscretization
    lₕ::Int
    Nₕ::Int
    basis::Basis
    mesh::OneDMesh{T}
    Rmax::T

    KohnShamSphericalDiscretization(lₕ::Int, basis::Basis, mesh::OneDMesh) = new{eltype(mesh)}(lₕ, length(basis), basis, mesh, last(mesh))
end

init_density_matrix(kd::KohnShamSphericalDiscretization)        = zeros(kd.lₕ+1, kd.Nₕ, kd.Nₕ), zeros(kd.lₕ+1, kd.Nₕ, kd.Nₕ)  
init_coeffs_discretization(kd::KohnShamSphericalDiscretization) = zeros(kd.lₕ+1, kd.Nₕ, kd.Nₕ)
init_energy(kd::KohnShamSphericalDiscretization)                = zeros(kd.lₕ+1, kd.Nₕ)
init_occupation(kd::KohnShamSphericalDiscretization)            = zeros(kd.lₕ+1, kd.Nₕ)

function build_kinetic!(kd::KohnShamSphericalDiscretization, Kin, A, M₋₂)
    @unpack lₕ = kd
    for l ∈ 0:lₕ
        @. Kin[l+1,:,:] = - A - l*(l+1)*M₋₂
    end 
end

function build_coulomb!(kd::KohnShamSphericalDiscretization, Coul, model, M₋₁)
    @unpack lₕ = kd
    for l ∈ 0:lₕ
        Coul[l+1,:,:] .= - charge(model) .* M₋₁
    end 
end

function build_hartree(kd::KohnShamSphericalDiscretization, Hartree)


end

function build_exchange_corr!(kd::KohnShamSphericalDiscretization, exc_mat, ρ, exc::AbstractExchangeCorrelation; quad_method, quad_retol, quad_atol)
    @unpack Nₕ, basis, Rmax = kd
    if !isthereExchangeCorrelation(exc)
        return exc_mat .= zeros(Nₕ,Nₕ)
    end
    for i ∈ eachindex(basis)
        for j ∈ eachindex(basis)
            if j<i
                exc_mat[i,j] = exc_mat[j,i]
            else
                exc_mat[i,j] = approximate_integral(x -> exc.vxc(ρ(x)), (0, Rmax) ; method = quad_method, retol = quad_retol, abstol = quad_atol)
            end
        end
    end
end

function build_density!(kd::KohnShamSphericalDiscretization, Dstar, U, n)
    @unpack lₕ, Nₕ, basis, mesh = kd
    density = zero_piecewiselaurantpolynomial(mesh)
    for l ∈ 1:lₕ+1
        for k ∈ 1:Nₕ
            if n[l,k] != 0
                eigen_vector = build_on_basis(basis, U[l,k,:])
                density += n[l,k] * eigen_vector * eigen_vector
            end
        end
    end
    density *= 1/4π
end

function build_energy()


end
