module KohnShamResolution

    # DEPENDANCES
    using LinearAlgebra
    using UnPack
    using Integrals
    using Memoize
    using HypergeometricFunctions
    using TensorOperations
    using Base.Threads
    using SparseArrays

    # ANNEXE
    include("utils.jl")
    include("computational tools.jl")

    # MESH
    export Mesh, linmesh, geometricmesh
    include("mesh.jl")

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

    
    # FEM
    abstract type Basis end
    abstract type AbstractLaurentPolynomialBasis <: Basis end 

    export LaurentPolynomialBasis
    export mass_matrix, stiffness_matrix, weight_mass_matrix, weight_mass_vector, vector_mass_matrix, vectorweight_mass_matrix,
           weight_mass_3tensor

    # SHORT BASIS
    include("fem/shortbasis/utils_computations.jl")

    #### OLD
    export AbstractShortElements
    export getpolynomial, getderivpolynomial
    include("fem/shortbasis/abstractshortelement.jl")

    #### NEW
    abstract type AbstractGenerator{T} end
    @inline Base.eltype(::AbstractGenerator{T}) where T = T
    @inline Base.length(gen::AbstractGenerator) = gen.size
    @inline getpolynomial(gen::AbstractGenerator, n::Int) = gen[n]
    @inline getderivpolynomial(gen::AbstractGenerator, n::Int) = getderivpolynomial(gen)[n]

    #### OLD
    export InfoElement
    export ShortPolynomialBasis
    export bottom_type
    include("fem/shortbasis/shortbasis.jl")

    #### NEW
    export PolynomialBasis
    include("fem/newbasis/basis.jl")
    include("fem/newbasis/matrices.jl")


    #### OLD
    export P1Elements, ShortP1Basis, IntLegendreElements, ShortIntLegendreBasis, DiffLegendreElements, ShortDiffLegendreBasis
    include("fem/shortbasis/elements.jl")

    #### NEW
    export IntLegendreGenerator, P1IntLegendreGenerator
    include("fem/newbasis/generators.jl")

    #### OLD

    include("fem/shortbasis/fill_mass_matrix.jl")

    export InfoBlock, CombineShortPolynomialBasis
    include("fem/shortbasis/combineshortbasis.jl")

    export ShortP1IntLegendreBasis
    include("fem/shortbasis/combined elements.jl")

    
    # KOHN-SHAM MODEL
    export ExchangeCorrelation,NoExchangeCorrelation, SlaterXα, LSDA
    export exc, vxc, vxcUP, vxcDOWN
    export KohnShamExtended, ReducedHartreeFock
    include("models.jl")
    
    # SOLVER &CO
    export LogConfig, LogBook
    include("log.jl")

    export Solver, SolverOptions
    include("solver.jl")

    abstract type KohnShamDiscretization end

    export LDADiscretization, LSDADiscretization
    include("lda_discretization.jl")
    include("lsda_discretization.jl")

    abstract type SCFMethod end

    export DFTProblem
    include("problem.jl")

    export aufbau!
    include("aufbau.jl")

    export CDA
    include("methods.jl")

    export KohnShamSolution, eigenvector, density
    include("solution.jl")

    export groundstate
    include("groundstate.jl")
end
