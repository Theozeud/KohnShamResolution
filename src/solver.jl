struct SolverOptions{T,Teps}
    lₕ::T   
    Nₕ::T    
    #Rmax::Tr    
    #Nr::T
    #Ngr::T
    #Ngθ::T
    ε::Teps      
end


struct KhonShamSolver{Tn,TVC}     
    cache::AbstractKohnShamCache
    opts::SolverOptions
    ϵ
    U
    n
    R
    Rprev
    niter::Tn
    val_crit::TVC
    current_crit
end