module KohnShamResolution

    # DEPENDANCIES
    using LinearAlgebra
    using SparseArrays
    using FillArrays
    using BlockDiagonals
    using LinearOperators
    using TensorOperations
    using Krylov
    using Optim
    using Integrals
    using HypergeometricFunctions
    using UnPack
    using Memoize
    using Base.Threads
    
    # ANNEXE
    include("utils.jl")
    include("maths.jl")

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
    abstract type KohnShamDiscretization end

    abstract type SCFMethod end
    abstract type SCFCache end
    abstract type SCFSolution end



    export LogConfig, LogBook
    include("log.jl")

    export KohnShamSolver, SolverOptions
    include("solver.jl")


    loopheader!(::SCFCache, ::SCFMethod, ::KohnShamSolver)      = nothing
    performstep!(::SCFCache, m::SCFMethod, ::KohnShamSolver)    = @warn "No performstep for the method $(typeof(m))"
    loopfooter!(::SCFCache, ::SCFMethod, ::KohnShamSolver)      = nothing
    monitor(::SCFCache, ::SCFMethod, ::KohnShamSolver)          = nothing
    register!(::SCFCache, ::SCFMethod, ::KohnShamSolver)        = nothing
    create_cache_method(m::SCFMethod, 
                        ::KohnShamDiscretization)               = @warn "No creation of cache for the method $(typeof(m))"
    switch!(c2::SCFCache, c1::SCFCache)                         = @warn "No way to init $(typeof(c2)) from of cache for the method ($(typeof(c1))"
    

    export LDADiscretization, LSDADiscretization
    include("discretization/lda.jl")
    include("discretization/lsda.jl")

    ###
    # TO REMOVE I GUESS ??? Not Sure!
    export DFTProblem
    include("problem.jl")
    ###
    
    export CDA, ODA, Quadratic
    include("methods/rca.jl")
    include("methods/oda_procedure.jl")
    #include("methods/quadratic.jl")
    include("methods/combined.jl")

    export aufbau!
    include("aufbau.jl")

    export KohnShamSolution, eigenvector, density, total_charge
    include("solution.jl")

    export groundstate
    include("groundstate.jl")
end
