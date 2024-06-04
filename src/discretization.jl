abstract type KohnShamDiscretization end


struct KohnShamSphericalDiscretization <: KohnShamDiscretization
    lₕ::Int
    Nₕ::Int
    basis::Basis
end

init_density_matrix(kd::KohnShamSphericalDiscretization) = zeros(kd.lₕ+1, kd.Nₕ, kd.Nₕ), zeros(kd.lₕ+1, kd.Nₕ, kd.Nₕ)  
init_coeffs_discretization(kd::KohnShamSphericalDiscretization) = zeros(kd.lₕ+1, (2*kd.lₕ+1)*kd.Nₕ, kd.Nₕ)
init_energy(kd::KohnShamSphericalDiscretization) = zeros(kd.lₕ+1,(2*kd.lₕ+1)*kd.Nₕ)
init_occupation(kd::KohnShamSphericalDiscretization) = zeros(kd.lₕ+1,(2*kd.lₕ+1)*kd.Nₕ)


function build_density_star!(kd::KohnShamSphericalDiscretization, temp_Dstar, temp_U, temp_n)
    @unpack lₕ, Nₕ = kd
    for l ∈ 1:lₕ+1
        for k ∈ 1:(2*(l-1)+1)*Nₕ
            temp_Dstar[l] .+= temp_n[l,k]*tensorproduct(temp_U[l,k], temp_U[l,k])
        end
    end
end

function build_kinetic!(kd::KohnShamSphericalDiscretization, Kin, A, M₋₂)
    @unpack lₕ = kd
    for l ∈ 0:lₕ
        Kin[l+1,:,:] .= - A .- l*(l+1)*M₋₂
    end 
end

function build_coulomb!(kd::KohnShamSphericalDiscretization, Coul, model, M₋₁)
    @unpack lₕ = kd
    for l ∈ 0:lₕ
        Coul[l+1,:,:] .= - 2*charge(model)*(2*l+1) .* M₋₁
    end 
end

function build_exchange_corr()


end

function build_potential()


end