struct SolverOptions{T}
    scftol::T       # SCF tolerance
    maxiter::Int    # Maximum of iteration done
    quad_method     # Method to use for quadrature of integrals
    quad_reltol::T  # Relativ tolerance for the quadrature of integrals
    quad_abstol::T  # Absolute tolerance for the quadrature of integrals
    hartree::T      # Coefficient multiply to the Hartree Matrix : 0 -> no hartree term, 1-> full hartree term
    degen_tol::T    
end

mutable struct KhonShamSolver
    niter                           # Number of iterations
    stop_criteria                   # Current stopping criteria
    discretization                  # Discretization parameters
    model                           # Model
    method                          # Iterative method
    opts                            # Solver options
    D                               # Density Matrix at current time
    Dprev                           # Density Matrix at previous time
    U                               # Coefficient of orbitals at current time
    ϵ                               # Orbitals energy at current time
    n                               # Occupation number at current time 
    energy                          # Total energy at current time
    logbook                         # LogBook
end