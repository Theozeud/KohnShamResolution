```
    __groundstate()
    
```
function groundstate(model::AbstractDFTModel, method::AbstractKohnShamResolutionMethod; kwargs...)
    solver = __init(model, method; kwargs...)
    __solve!(solver, method)
    __makesolution(solver)
end

```
    __init()
    
```
function __init(model::AbstractDFTModel, method::AbstractKohnShamResolutionMethod; 
    mesh::AbstractMesh, ci = 0.0)

    # All Checks
    if !ismethod_for_model(method, model)
        error("This method can't be used for this model")
    end

    # Init Cache
    cache = init_cache(method)

    solver(mesh, cache, ci)
end

```
    __solve!()
    
```
function __solve!(solver::KhonShamSolver, method::AbstractKohnShamResolutionMethod)
    while stopping_criteria(method, solver)
        #loopheader(solver)
        performstep!(method, solver)
        #loopfooter(solver)
    end
end

```
    __makesolution()
    
```
function __makesolution(solver::KhonShamSolver)

end
