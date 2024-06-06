struct SolverOptions{Teps}
    ε::Teps 
    maxiter::Int     
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