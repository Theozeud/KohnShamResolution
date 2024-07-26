
# Solve the eigen problem
function _eigen(basis, ρ, Rmin, Rmax; l = 0, ϵ = 1, z = 1)
    deriv_basis = deriv(basis)
    A   = Symmetric(mass_matrix(deriv_basis))
    M₀  = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(weight_mass_matrix(basis, -1))
    # Compute Hartree matrix
    F =  weight_mass_vector(basis, x->4π*ρ(x) * x)
    Cᵨ = 4π * integrate(ρ * Monomial(2), Rmin, Rmax)
    Hart = vectorweight_mass_matrix(basis, A\F, Monomial(-1))+ Cᵨ/(Rmax-Rmin) * M₀
    if l == 0
        H = T(0.5) * A  -  z .* M₋₁ + ϵ * Hart
        return  eigen(H, M₀)
    else
        M₋₂ = Symmetric(weight_mass_matrix(basis, -2))
        H = T(0.5) * A  -  z.* M₋₁ + T(0.5) * l*(l+1) * M₋₂ + ϵ * Hart
        return  eigen(H, M₀)
    end
end


function _build_density(basis, vect, z = 1)
    tmp = build_on_basis(basis, vect)
    1/4π * tmp*tmp*Monomial(-2) * z
end


function _relax(ρ, ρprev, t)
    ρ*t + (1-t)*ρprev
end

function _iter(basis, Rmin, Rmax; l = 0, ϵ = 1, tol = 1e-2, maxiter = 10, oda = 0.8, z = 1)
    ρprev = Monomial(0,0) 
    crit = 10*tol
    iter = 0
    arraycrit = []
    arrayvals = []
    while crit > tol && iter < maxiter
        @time "$crit" vals, vects = _eigen(basis, ρprev, Rmin, Rmax; l = l, ϵ = ϵ, z= z)
        ρ = _relax(_build_density(basis, vects[:,1], z), ρprev, oda)
        crit = _crit(ρ, ρprev, Rmin, Rmax)
        push!(arraycrit, crit)
        push!(arrayvals, min(vals...))
        ρprev = ρ
        iter += 1 
    end
    ρprev, arrayvals, arraycrit
end


function _crit(ρ, ρprev, Rmin, Rmax)
    X = range(Rmin, Rmax, 1000)
    sqrt((Rmax-Rmin)*sum(abs.(ρ.(X[1:end-1])-ρprev.(X[1:end-1])).^2)/(length(X)-1))
end

