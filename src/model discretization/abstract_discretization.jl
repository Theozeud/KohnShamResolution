abstract type KohnShamDiscretization end

init_density_matrix(::KohnShamDiscretization) = @error "Missing functions"
init_coeffs_discretization(::KohnShamDiscretization) = @error "Missing functions"
init_energy(::KohnShamDiscretization) = @error "Missing functions"
init_occupation(::KohnShamDiscretization) = @error "Missing functions"

build_density_star!(::KohnShamDiscretization, args...; kwargs...) = @error "Missing functions"
build_kinetic!(::KohnShamDiscretization, args...; kwargs...) = @error "Missing functions"
build_coulomb!(::KohnShamDiscretization, args...; kwargs...) = @error "Missing functions"
build_exchange_corr(::KohnShamDiscretization, args...; kwargs...) = @error "Missing functions"
build_potential(::KohnShamDiscretization, args...; kwargs...) = @error "Missing functions"


abstract type AbstractKohnShamCache end