module KohnShamResolution

    # DEPENDANCES
    using LinearAlgebra
    using UnPack
    using Integrals
    using Memoize
    using HypergeometricFunctions
    using TensorOperations
    using Base.Threads

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

    
    # FINITE ELEMENT BASIS
    abstract type Basis end
    abstract type AbstractLaurentPolynomialBasis <: Basis end 


    export LaurentPolynomialBasis
    export mass_matrix, stiffness_matrix, weight_mass_matrix, weight_mass_vector, vector_mass_matrix, vectorweight_mass_matrix,
           weight_mass_3tensor
    export build_on_basis

    # SHORT BASIS
    include("finite element basis/shortbasis/utils_computations.jl")

    export AbstractShortElements
    export getpolynomial, getderivpolynomial
    include("finite element basis/shortbasis/abstractshortelement.jl")

    #### NEW
    abstract type AbstractGenerator{T} end
    @inline Base.eltype(::AbstractGenerator{T}) where T = T
    @inline Base.length(gen::AbstractGenerator) = gen.size
    @inline getpolynomial(gen::AbstractGenerator, n::Int) = gen[n]
    @inline getderivpolynomial(gen::AbstractGenerator, n::Int) = getderivpolynomial(gen)[n]

    export InfoElement
    export ShortPolynomialBasis, build_basis
    export bottom_type
    include("finite element basis/shortbasis/shortbasis.jl")

    #### NEW
    export PolynomialBasis
    include("finite element basis/newbasis/basis.jl")

    export P1Elements, ShortP1Basis, IntLegendreElements, ShortIntLegendreBasis, DiffLegendreElements, ShortDiffLegendreBasis
    include("finite element basis/shortbasis/elements.jl")

    #### NEW
    export IntLegendreGenerator
    include("finite element basis/newbasis/generators.jl")

    include("finite element basis/shortbasis/fill_mass_matrix.jl")

    export InfoBlock, CombineShortPolynomialBasis
    include("finite element basis/shortbasis/combineshortbasis.jl")

    export ShortP1IntLegendreBasis
    include("finite element basis/shortbasis/combined elements.jl")

    
    # KOHN-SHAM MODEL
    export AbstractExchangeCorrelation, ExchangeCorrelation, NoExchangeCorrelation, KohnShamExtended
    export build_SlaterXα, exc_SlaterXα, vxc_SlaterXα
    export ReducedHartreeFock, SlaterXα
    include("models.jl")
    
    # SOLVER &CO
    export LogConfig, LogBook
    include("log.jl")

    include("solver.jl")

    export KohnShamRadialDiscretization
    include("radial_discretization.jl")

    abstract type SCFMethod end

    export DFTProblem
    include("problem.jl")

    export aufbau!
    include("aufbau.jl")

    export CDA
    include("methods.jl")

    export KohnShamSolution
    include("solution.jl")

    export groundstate
    include("groundstate.jl")
end
