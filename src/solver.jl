struct SolverOptions{T}
    ε::T
    maxiter::Int
    quad_method
    quad_reltol::T
    quad_abstol::T
    hartree::T
    degen_tol::T
    potential::Symbol
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
    Ehisto
    ϵhisto
    Energyhisto
    niter::Int
    values_stop_crit::Vector{TSC}
    current_stop_crit::TSC
    opts::SolverOptions
end