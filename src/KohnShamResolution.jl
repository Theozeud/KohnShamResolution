module KohnShamResolution

    using LinearAlgebra
    using SparseArrays
    using TimerOutputs
    
    
    include("utils.jl")

    include("computation_tools.jl")

    export OneDMesh
    export mesh
    
    include("mesh.jl")

    export KohnShamExtended
    export exchcorr, charge, nbelec, potential

    include("model.jl")

    include("abstractmethods.jl")

    include("solver.jl")

    export ODA

    include("methods.jl")


    include("solution.jl")

    export groundstate

    include("groundstate.jl")
    

end
