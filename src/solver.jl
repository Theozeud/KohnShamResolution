struct SolverOptions{T,Teps}
    ε::Teps      
end


struct KhonShamSolver{TSC}
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