struct SolverOptions{T,Teps}
    lₕ::T   
    Nₕ::T    
    #Rmax::Tr    
    #Nr::T
    #Ngr::T
    #Ngθ::T
    ε::Teps      
end


struct KhonShamSolver{TSC}
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