module KohnShamResolution

    using LinearAlgebra
    using SparseArrays
    using TimerOutputs
    using ProgressMeter
    using UnPack
    using Integrals
    
    include("utils.jl")
    include("tools/computation_tools.jl")

    export LaurentPolynomial
    export Monomial, deg, degmax, degmin, haslog, ismonomial, iszero
    export integrate!, integrate, deriv!, deriv, scalar_product, elag!
    include("basis/laurentpolynomial.jl")

    export LaurentPolynomialBasis
    export mass_matrix, weight_mass_matrix, build_on_basis
    include("basis/laurentpolynomialbasis.jl")

    export OneDMesh
    export mesh, find_index
    include("mesh.jl")

    export PiecewiseLaurentPolynomial
    export get_support
    include("basis/piecewiselaurentpolynomial.jl")

    export HatFunctionP1, HatBasis
    include("basis/hat_functions.jl")

    export AbstractExchangeCorrelation, ExchangeCorrelation, NoExchangeCorrelation, KohnShamExtended
    export build_SlaterXα, exc_SlaterXα, vxc_SlaterXα
    export exchcorr, charge, nbelec, potential
    export ReducedHartreeFock, SlaterXα
    include("model/model.jl")

    export KohnShamDiscretization, KohnShamSphericalDiscretization
    include("discretization/abstract_discretization.jl")
    include("discretization/spherical_discretization.jl")

    include("methods/abstractmethods.jl")

    export DFTProblem
    include("problem.jl")

    include("solver.jl")

    export ODA
    include("methods/oda.jl")

    include("solution.jl")

    export groundstate
    include("groundstate.jl")
    

end
