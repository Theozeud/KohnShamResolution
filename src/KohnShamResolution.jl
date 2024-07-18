module KohnShamResolution

    using LinearAlgebra
    using UnPack
    using Integrals
    using Memoize
    using LambertW
    using HypergeometricFunctions

    using SparseArrays
    using TimerOutputs
    using ProgressMeter

    ########
    # ANNEXE
    include("utils.jl")
    include("computational tools.jl")

    ######
    # MESH
    export OneDMesh
    export mesh, find_index, linmesh, LogRange, logmesh, LinearExpRange, linearexpmesh
    include("mesh.jl")

    ####################
    # LAURENT POLYNOMIAL
    export AbstractPolynomial
    include("laurentpolynomial/abstract polynomial.jl")

    export LaurentPolynomial
    export Polynomial, Monomial, RandMonomial, RandPolynomial, RootsPolynomial
    export deg, degmax, degmin, haslog, ismonomial, iszero
    export integrate!, integrate, deriv!, deriv, scalar_product, normL2, elag!, diveucl
    include("laurentpolynomial/laurentpolynomial.jl")

    include("laurentpolynomial/integration_tools.jl")

    export RationalFraction, SommeRationalFraction
    export fraction_decomp, inEntpart
    include("laurentpolynomial/rationalfraction.jl")

    export Legendre, intLegendre
    include("laurentpolynomial/legendre polynomial.jl")

    export PiecewiseLaurentPolynomial
    export get_support
    include("laurentpolynomial/piecewiselaurentpolynomial.jl")

    ######################
    # FINITE ELEMENT BASIS
    export Basis, AbstractLaurentPolynomialBasis
    include("finite element basis/abstract polynomial basis.jl")
    
    export LaurentPolynomialBasis
    export mass_matrix, weight_mass_matrix, build_on_basis
    include("finite element basis/completebasis/laurentpolynomialbasis.jl")

    export HatBasis, P1Basis, BubbleBasis, IntLegendreBasis
    export HatFunctionP1, FunctionP2_node, FunctionP2_mid, P2Basis, IntLegendre_element
    include("finite element basis/completebasis/complete basis.jl")

    # SHORT BASIS
    include("finite element basis/shortbasis/utils_computations.jl")

    export AbstractShortElements
    export getpolynomial, isnormalized
    include("finite element basis/shortbasis/abstractshortelement.jl")

    export InfoElement
    export ShortPolynomialBasis, build_basis
    export bottom_type
    include("finite element basis/shortbasis/shortbasis.jl")

    export DefaultElements, P1Elements, ShortP1Basis, IntLegendreElements, ShortIntLegendreBasis, DiffLegendreElements, ShortDiffLegendreBasis
    include("finite element basis/shortbasis/elements.jl")

    include("finite element basis/shortbasis/fill_mass_matrix.jl")

    export InfoBlock, CombineShortPolynomialBasis
    include("finite element basis/shortbasis/combineshortbasis.jl")

    export ShortP1IntLegendreBasis
    include("finite element basis/shortbasis/combine elements.jl")

    #################
    # KOHN-SHAM MODEL
    export AbstractExchangeCorrelation, ExchangeCorrelation, NoExchangeCorrelation, KohnShamExtended
    export build_SlaterXα, exc_SlaterXα, vxc_SlaterXα
    export exchcorr, charge, nbelec, potential
    export ReducedHartreeFock, SlaterXα
    include("model/kohnsham_model.jl")

    ##########################
    # KOHN-SHAM DISCRETIZATION
    export KohnShamDiscretization, KohnShamSphericalDiscretization
    include("model discretization/abstract_discretization.jl")
    include("model discretization/spherical_discretization.jl")

    ###################
    # KOHN-SHAM METHODS
    include("methods/abstractmethods.jl")

    ###############
    # SOLVER AND CO
    export DFTProblem
    include("problem.jl")

    include("solver.jl")

    export ODA, ConstantODA
    include("methods/oda.jl")

    export KohnShamSolution
    export eigenvalue, eigenvector, occup
    include("solution.jl")

    export groundstate
    include("groundstate.jl")
end
