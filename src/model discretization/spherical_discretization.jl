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


function build_kinetic!(::KohnShamSphericalDiscretization, Kin, A)
    @. Kin[1,:,:] =  1/2 * A
end

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

function build_hartree_deprecated!(kd::KohnShamSphericalDiscretization, Hartree, ρ)
    @unpack basis, Rmin, Rmax = kd
    int1 = integrate(Monomial(1) * ρ)
    int2 = integrate(Monomial(2) * ρ)
    potential   = 4π*(int1(Rmax) - int1 + Monomial(-1) * (int2 - int2(Rmin)))
    ∇potential  = - Monomial(-2) * (int2 - int2(Rmin))
    for i ∈ eachindex(basis)
        for j ∈ eachindex(basis)
            if j<i
                Hartree[i,j] = Hartree[j,i]
            else
                #Hartree[i,j] = integrate(potential * basis[i] * basis[j], Rmin, Rmax)
                #Hartree[i,j] = approximate_integral(x -> potential(x) * basis[i](x) * basis[j](x), (Rmin, Rmax) ; method = QuadGKJL(), reltol = 1e-3, abstol = 1e-3)
                #int_ᵢⱼ = integrate(basis[i] * basis[j])
                int_ᵢⱼ = integrate(build_basis(basis, i) * build_basis(basis, j))
                M₀ᵢⱼ = int_ᵢⱼ - int_ᵢⱼ(Rmin)
                Hartree[i,j] = potential(Rmax)*M₀ᵢⱼ(Rmax) - integrate(∇potential * M₀ᵢⱼ, Rmin, Rmax)
            end
        end
    end
    Hartree
end

function build_hartree!(kd::KohnShamSphericalDiscretization, Hartree, ρ, opt)
    @unpack basis, Rmin, Rmax = kd
    if opt == :integral
        build_hartree_integral!(kd::KohnShamSphericalDiscretization, Hartree, ρ)
    elseif opt == :pde
        deriv_basis = deriv(basis)
        A   = mass_matrix(deriv_basis)
        M₀ = mass_matrix(basis)
        build_hartree_pde!(kd, Hartree, ρ, A, M₀)
    end
    Hartree
end

function build_hartree_integral!(kd::KohnShamSphericalDiscretization, Hartree, ρ)
    @unpack basis, Rmin, Rmax = kd
    potential = nothing
    fill_weight_mass_matrix!(basis, potential, Hartree)
    Hartree
end

function build_hartree_pde!(kd::KohnShamSphericalDiscretization, Hartree, ρ, A, M₀)
    @unpack basis, Rmin, Rmax = kd
    rho(x) = ρ(x)
    f(x) = 4π * rho(x) * x
    F =  weight_mass_vector(basis, f)
    coeff = A\F
    Cᵨ = 4π * integrate(ρ * Monomial(2), Rmin, Rmax)
    MV  = vectorweight_mass_matrix(basis, coeff, Monomial(-1)) 
    @. Hartree = MV + Cᵨ/(Rmax-Rmin) * M₀
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

function build_eigenvector(kd::KohnShamSphericalDiscretization, U; Index = CartesianIndices((1:lₕ+1, 1:Nₕ)))
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
