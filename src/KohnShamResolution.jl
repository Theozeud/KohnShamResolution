module KohnShamResolution

    using LinearAlgebra
    using SparseArrays
    using UnPack
    using Integrals
    using Memoize

    using TimerOutputs
    using ProgressMeter
    
    include("utils.jl")
    include("computational tools.jl")

    export LaurentPolynomial
    export Polynomial, Monomial, deg, degmax, degmin, haslog, ismonomial, iszero
    export integrate!, integrate, deriv!, deriv, scalar_product, elag!
    include("laurentpolynomial/laurentpolynomial.jl")

    export Legendre
    include("laurentpolynomial/legendre polynomial.jl")

    export OneDMesh
    export mesh, find_index, linmesh, LogRange, logmesh
    include("mesh.jl")

    export PiecewiseLaurentPolynomial
    export get_support
    include("laurentpolynomial/piecewiselaurentpolynomial.jl")

    export LaurentPolynomialBasis
    export mass_matrix, weight_mass_matrix, build_on_basis
    include("finite element basis/laurentpolynomialbasis.jl")

    export HatBasis, P1Basis, BubbleBasis, IntLegendreBasis
    export HatFunctionP1, FunctionP2_node, FunctionP2_mid, P2Basis, IntLegendre_element
    include("finite element basis/usual basis.jl")

    export AbstractExchangeCorrelation, ExchangeCorrelation, NoExchangeCorrelation, KohnShamExtended
    export build_SlaterXα, exc_SlaterXα, vxc_SlaterXα
    export exchcorr, charge, nbelec, potential
    export ReducedHartreeFock, SlaterXα
    include("model/kohnsham_model.jl")

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
