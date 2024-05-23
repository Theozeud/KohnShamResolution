#=
    __groundstate()
    
=#
function groundstate(model::AbstractDFTModel, method::AbstractKohnShamResolutionMethod; kwargs...)
    solver = __init(model, method; kwargs...)
    __solve!(solver, method)
    __makesolution(solver)
end


#=
    __init()
    
=#
function __init(model::AbstractDFTModel, method::AbstractKohnShamResolutionMethod; 
    mesh::AbstractMesh, ci = 0.0, lₕ, Nₕ, _ϵ)

    # All Checks
    if !ismethod_for_model(method, model)
        error("This method can't be used for this model")
    end

    # Init Cache
    cache = init_cache(method)

    # Registering SolverOptions
    opts = SolverOptions(lₕ,Nₕ,_ϵ)

    # Init storage array
    ϵ = spzeros(lₕ+1,(2lₕ+1)Nₕ)
    U = spzeros(lₕ+1, (2lₕ+1)Nₕ, Nₕ)
    n = spzeros(lₕ+1, (2lₕ+1)Nₕ)
    R = spzeros(lₕ+1, Nₕ, Nₕ)
    Rprev = spzeros(lₕ+1, Nₕ, Nₕ)

    solver(mesh, cache, opts, ϵ, U, n, R, Rprev)
end

#=
    __solve!()
    
=#
function __solve!(solver::KhonShamSolver, method::AbstractKohnShamResolutionMethod)
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
    solver.Rprev = solver.R
end



#=
    __makesolution()
    
=#
function __makesolution(solver::KhonShamSolver)

end
