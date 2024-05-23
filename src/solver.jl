struct SolverOptions{T,Tr,Teps}
    lₕ::T   
    Nₕ::T    
    #Rmax::Tr    
    #Nr::T
    #Ngr::T
    #Ngθ::T
    ε::Teps      
end


struct KhonShamSolver{TS}
    mesh::AbstractMesh          
    cache::AbstractKohnShamCache
    opts::SolverOptions
    ϵ::AbstractMatrix
    U::AbstractMatrix
    n::AbstractMatrix
    R::AbstractMatrix
    Rprev::AbstractMatrix
end