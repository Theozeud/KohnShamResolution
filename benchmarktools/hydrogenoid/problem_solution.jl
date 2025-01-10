# STRUCTURE HYDROGENOID PROBLEM

struct HydrogenoidOption
    nλ         # Number of the eigenvalue wanted
    nU         # Number of the eigenvector wanted
end

struct HydrogenoidProblem
    T           # Type of numbers
    z           # Number of electrons
    l           # 2nd quantum number of orbital
    Rmax        # Spatial cut-off
    Nmesh       # Number of points of the discretization
    typemesh    # Type of the Mesh
    optsmesh    # Options for the Mesh
    typebasis   # Type of the basis
    optsbasis   # Options for the basis
    opts        # Option for for the solution
    name        # Name of the problem

    function HydrogenoidProblem(; 
        T, z, l, Rmax, Nmesh,typemesh, typebasis, optsmesh, optsbasis, 
        nλ = 1:10, nU = nλ, name = "")
        _nU = isnothing(nU) ? (1:0) : nU
        opts = HydrogenoidOption(nλ, _nU)
        new(T, z, l, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis, opts, name)
    end

    function HydrogenoidProblem(prob; 
        T = prob.T, z = prob.z, l = prob.l, Rmax = prob.Rmax, Nmesh = prob.Nmesh,
        typemesh = prob.typemesh, typebasis = prob.typebasis, 
        optsmesh = prob.optsmesh, optsbasis = prob.optsbasis, 
        nλ = prob.nλ, nU = prob.nU, name = prob.name)
        _nU = isnothing(nU) ? (1:0) : nU
        opts = HydrogenoidOption(nλ, _nU)
        new(T, z, l, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis, opts, name)
    end
end



# STRUCTURE HYDROGENOID SOLUTION

struct HydrogenoidSolution
    problem         # Original problem
    nλ              # Number of the eigenvalue computed
    nU              # Number of the eigenvector computed
    λ               # Numerical eigenvalue
    U               # Numerical eigenvector
    λtheo           # Theoretical eigenvalue
    Utheo           # Theoretical eigenvector
    Δλ              # |λ - λtheo|
    ΔU              # norm(U - Utheo) TO CLARIFY THE COMPUTATION OF THE NORM

    function HydrogenoidSolution(problem, λ, U)
        nλ = problem.opts.nλ
        nU = problem.opts.nU
        _λ = λ[nλ]
        _U = nothing
        λtheo = theoretical_eigenvalue(problem)
        Utheo =  theoretical_eigenvector(problem)
        Δλ = abs.(_λ .- λtheo)
        ΔU = Δλ      ## TO DO
        new(problem, nλ, nU, _λ, _U, λtheo, Utheo, Δλ, ΔU)
    end
end 

# STRUCTURE HYDROGENOID FOR ANALYSING CONVERGENCE

struct HydrogenoidConvergenceNmesh
    probs           # Set of problems
    vecNmesh        # Set of Nmesh used
    ΔΛ              # Dict of errors on eigenvalues : for each problem,
                    # there is a vector of errors depending on Nmesh
    ΔU              # Dict of errors on eigenvectors : for each problem,
                    # there is a vector of errors depending on Nmesh
    num             # Number of eigenvalue and eigenvector used 
end

struct HydrogenoidConvergenceRmax
    probs           # Set of problems
    vecRmax         # Set of Rmax used
    ΔΛ              # Dict of errors on eigenvalues : for each problem,
                    # there is a vector of errors depending on Rmax
    ΔU              # Dict of errors on eigenvectors : for each problem,
                    # there is a vector of errors depending on Rmax
    num             # Number of eigenvalue and eigenvector used
end

