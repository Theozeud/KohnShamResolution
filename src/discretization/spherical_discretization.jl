struct KohnShamSphericalDiscretization <: KohnShamDiscretization
    lₕ::Int
    Nₕ::Int
    basis::Basis
    mesh::OneDMesh

    KohnShamSphericalDiscretization(lₕ::Int, basis::Basis, mesh::OneDMesh) = new(lₕ, length(basis), basis, mesh)
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

function build_exchange_corr()


end

function build_density!(kd::KohnShamSphericalDiscretization, temp_Dstar, temp_U, temp_n)
    @unpack lₕ, Nₕ, basis, mesh = kd
    density = zero_piecewiselaurantpolynomial(mesh)
    for l ∈ 1:lₕ+1
        for k ∈ 1:Nₕ
            if temp_n[l,k] != 0
                eigen_vector = build_on_basis(basis, temp_U[l,k,:])
                density += temp_n[l,k] * eigen_vector * eigen_vector
            end
        end
    end
    density *= 1/4π
end

