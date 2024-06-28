struct SolverOptions{T}
    ε::T
    maxiter::Int
    quad_method
    quad_reltol::T
    quad_abstol::T
    hartree::T
    degen_tol::T
end

mutable struct KhonShamSolver{TSC}
    discretization
    model
    D
    Dprev
    U  
    ϵ
    n
    niter::Int
    values_stop_crit::Vector{TSC}
    current_stop_crit::TSC
    cache::AbstractKohnShamCache
    opts::SolverOptions
end