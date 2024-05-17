struct SolverOptions{T,Tr,Teps}
    lh::T      # Cut-off sections of the decomposition in eigenspace
    Rmax::Tr    # Cut-off of the radial domain
    Nr::T
    Ngr::T
    Ngθ::T
    ϵ::Teps      
end


struct KhonShamSolver{TS}
    mesh::AbstractMesh
    cache::AbstractKohnShamCache
    sol::TS
    opts::SolverOptions

end