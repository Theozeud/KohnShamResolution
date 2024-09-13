struct SolverOptions{T}
    ε::T            # Alg tolerance
    maxiter::Int    # Maximum of iteration done
    quad_method     # Method to use for quadrature of integrals
    quad_reltol::T  # Relativ tolerance for the quadrature of integrals
    quad_abstol::T  # Absolute tolerance for the quadrature of integrals
    hartree::T      # coefficient multiply to the Hartree Matrix : 0 -> no hartree term, 1-> full hartree term
    degen_tol::T    
    light::Bool
end

mutable struct KhonShamSolver{TSC}
    discretization
    model
    method
    D
    Dprev
    U  
    ϵ
    n
    ϵhisto
    Energyhisto
    niter::Int
    values_stop_crit::Vector{TSC}
    current_stop_crit::TSC
    opts::SolverOptions
end