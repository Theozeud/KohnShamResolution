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
