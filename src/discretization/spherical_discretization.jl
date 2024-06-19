struct KohnShamSphericalDiscretization{T} <: KohnShamDiscretization
    lₕ::Int
    Nₕ::Int
    basis::Basis
    mesh::OneDMesh{T}
    Rmin::T
    Rmax::T

    KohnShamSphericalDiscretization(lₕ::Int, basis::Basis, mesh::OneDMesh) = new{eltype(mesh)}(lₕ, length(basis), basis, mesh, first(mesh), last(mesh))
end

init_density_matrix(kd::KohnShamSphericalDiscretization, T::Type)        = zero_piecewiselaurantpolynomial(kd.mesh, T), zero_piecewiselaurantpolynomial(kd.mesh, T) 
init_coeffs_discretization(kd::KohnShamSphericalDiscretization, T::Type) = zeros(T, kd.lₕ+1, kd.Nₕ, kd.Nₕ)
init_energy(kd::KohnShamSphericalDiscretization, T::Type)                = zeros(T, kd.lₕ+1, kd.Nₕ)
init_occupation(kd::KohnShamSphericalDiscretization, T::Type)            = zeros(T, kd.lₕ+1, kd.Nₕ)

function build_kinetic!(kd::KohnShamSphericalDiscretization, Kin, A, M₋₂)
    @unpack lₕ = kd
    for l ∈ 0:lₕ
        @. Kin[l+1,:,:] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
end

function build_coulomb!(kd::KohnShamSphericalDiscretization, Coul, model, M₋₁)
    @unpack lₕ = kd
    for l ∈ 0:lₕ
        Coul[l+1,:,:] .= - charge(model) .* M₋₁
    end 
end

function build_hartree!(kd::KohnShamSphericalDiscretization, Hartree, ρ)
    @unpack basis, Rmin, Rmax = kd
    int1 = integrate(Monomial(1) * ρ)
    int2 = integrate(Monomial(2) * ρ)
    potential = 4π*(int1 - int1(Rmin) + Monomial(-1) * (int2(Rmax) - int2))
    for i ∈ eachindex(basis)
        for j ∈ eachindex(basis)
            if j<i
                Hartree[i,j] = Hartree[j,i]
            else
                #Hartree[i,j] = integrate(potential * basis[i] * basis[j], Rmin, Rmax)
                Hartree[i,j] = approximate_integral(x -> potential(x) * basis[i](x) * basis[j](x), (Rmin, Rmax) ; method = QuadGKJL(), reltol = 1e-3, abstol = 1e-3)
            end
        end
    end
    Hartree
end

function build_exchange_corr!(kd::KohnShamSphericalDiscretization, exc_mat, ρ, exc::AbstractExchangeCorrelation; quad_method, quad_reltol, quad_abstol)
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

function build_eigenvector(kd::KohnShamSphericalDiscretization, U)
    @unpack lₕ, Nₕ, basis, mesh = kd
    eigen_vector = []
    for l ∈ 1:lₕ+1
        for k ∈ 1:Nₕ
            eiglk = build_on_basis(basis, U[l,:,k]) 
            push!(eigen_vector,  1/sqrt(4π)* eiglk / normL2(eiglk, mesh) * Monomial(-1)) 
        end
    end
    reshape(eigen_vector, (Nₕ, lₕ+1))
end


function build_density!(kd::KohnShamSphericalDiscretization, Dstar, U, n)
    @unpack lₕ, Nₕ, basis, mesh = kd
    for l ∈ 1:lₕ+1
        for k ∈ 1:Nₕ
            if n[l,k] != 0
                eigen_vector = build_on_basis(basis, U[l,:,k]) 
                eigen_vector = eigen_vector / normL2(eigen_vector, mesh) * Monomial(-1)
                Dstar += n[l,k] * eigen_vector * eigen_vector
            end
        end
    end
    Dstar *= 1/(4π)
end

function build_energy()


end
