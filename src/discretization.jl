abstract type KohnShamDiscretization end


struct KohnShamSphericalDiscretization <: KohnShamDiscretization
    lₕ::Int
    Nₕ::
    basis::Basis
end

init_density_matrix(::KohnShamSphericalDiscretization) = zeros(lₕ+1, Nₕ, Nₕ), zeros(lₕ+1, Nₕ, Nₕ)  
init_coeffs_discretization(::KohnShamSphericalDiscretization) = zeros(lₕ+1, (2lₕ+1)Nₕ, Nₕ)
init_energy(::KohnShamSphericalDiscretization) = zeros(lₕ+1,(2lₕ+1)Nₕ)
init_occupation(::KohnShamSphericalDiscretization) = zeros(lₕ+1,(2lₕ+1)Nₕ)


function build_density_star!(::KohnShamSphericalDiscretization, temp_Dstar, temp_U, temp_n)
    for l ∈ 1:lₕ+1
        for k ∈ 1:(2*(l-1)+1)*Nₕ
            temp_Dstar[l] .+= temp_n[l,k]*tensorproduct(temp_U[l,k], temp_U[l,k])
        end
    end
end
