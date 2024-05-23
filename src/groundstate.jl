#=
    groundstate()
    
=#
function groundstate(model::AbstractDFTModel, method::AbstractKohnShamResolutionMethod; kwargs...)
    solver = init(model, method; kwargs...)
    solve!(solver, method)
    makesolution(solver)
end


#=
    init()
    
=#
function init(model::AbstractDFTModel, method::AbstractKohnShamResolutionMethod; 
    mesh::AbstractMesh, ci = 0.0, lₕ, Nₕ, ε)

    # All Checks
    if !ismethod_for_model(method, model)
        error("This method can't be used for this model")
    end

    # Init Cache
    cache = init_cache(method)

    # Registering SolverOptions
    opts = SolverOptions(lₕ,Nₕ,ε)

    # Init storage array
    ϵ       = spzeros(lₕ+1,(2lₕ+1)Nₕ)
    U       = spzeros(lₕ+1, (2lₕ+1)Nₕ, Nₕ)
    n       = spzeros(lₕ+1, (2lₕ+1)Nₕ)
    R       = spzeros(lₕ+1, Nₕ, Nₕ)
    Rprev   = spzeros(lₕ+1, Nₕ, Nₕ)

    solver(mesh, cache, opts, ϵ, U, n, R, Rprev)
end

#=
    solve!()
    
=#
function solve!(solver::KhonShamSolver, method::AbstractKohnShamResolutionMethod)
    while stopping_criteria(method, solver, solver.opts.ϵ)
        loopheader!(solver)
        performstep!(method, solver.cache, solver)
        loopfooter!(solver)
    end
end


#=
    loopheader!()
    
=#
function loopheader!(::KhonShamSolver)
    
end


#=
    loopfooter!()
    
=#
function loopfooter!(solver::KhonShamSolver)
    solver.Rprev .= solver.R
end



#=
    __makesolution()
    
=#
function makesolution(solver::KhonShamSolver)

end
