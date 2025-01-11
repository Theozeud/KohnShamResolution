# Integrated Legendre Generators

struct IntLegendreGenerator{T} <: AbstractGenerator{T}
    polynomials::Vector{LaurentPolynomial{T}}
    derivpolynomials::Vector{LaurentPolynomial{T}}
    size::Int
    ordermin::Int
    ordermax::Int
    binf::T
    bsup::T

    function IntLegendreGenerator(T::Type = Float64; ordermin::Int = 2, ordermax = 2, binf::Real = -T(1), bsup::Real = T(1))
        @assert ordermin ≥ 1
        polynomials = LaurentPolynomial{T}[]
        derivpolynomials = LaurentPolynomial{T}[]  
        for n ∈ ordermin:ordermax
            Pₙ = Legendre(n-1; T = T, a = T(binf), b = T(bsup))
            push!(derivpolynomials, Pₙ)
            Qₙ = intLegendre(n-1; T = T, a = T(binf), b = T(bsup))
            push!(polynomials, Qₙ)
        end
        new{T}(polynomials, derivpolynomials, ordermax - ordermin + 1, ordermin, ordermax, T(binf), T(bsup))
    end
end

@inline Base.firstindex(::IntLegendreGenerator) = 1
@inline Base.eachindex(ilg::IntLegendreGenerator) = eachindex(ilg.polynomials)
@inline Base.getindex(ilg::IntLegendreGenerator, n::Int) =  ilg.polynomials[n] 
@inline getpolynomial(ilg::IntLegendreGenerator) = ilg.polynomials
@inline getderivpolynomial(ilg::IntLegendreGenerator) = ilg.derivpolynomials

function IntLegendreGenerator(mesh::Mesh, T::Type = Float64; kwargs...)
    generators = IntLegendreGenerator(T; kwargs...)
    size = generators.size * (lastindex(mesh) - 1)
    indices_cells       = zeros(Int, size, 1)
    indices_generators  = zeros(Int, size, 1)
    normalisation       = ones(Int, size)
    for i ∈ eachindex(mesh)[1:end-1]
        for n ∈ 1:generators.size
            indices_cells[(i-1) * generators.size + n, 1]       = i
            indices_generators[(i-1) * generators.size + n, 1]  = n
        end
    end
    PolynomialBasis(generators, mesh, size, indices_cells, indices_generators, normalisation) 
end


# Integrated Legendre Generators

struct P1IntLegendreGenerator{T} <: AbstractGenerator{T}
    polynomials::Vector{LaurentPolynomial{T}}
    derivpolynomials::Vector{LaurentPolynomial{T}}
    size::Int
    ordermin::Int
    ordermax::Int
    binf::T
    bsup::T

    function P1IntLegendreGenerator(T::Type = Float64; ordermin::Int = 2, ordermax = 2, binf::Real = -T(1), bsup::Real = T(1))
        @assert ordermin ≥ 1
        polynomials = LaurentPolynomial{T}[]
        derivpolynomials = LaurentPolynomial{T}[]  
        for n ∈ ordermin:ordermax
            Pₙ = Legendre(n-1; T = T, a = T(binf), b = T(bsup))
            push!(derivpolynomials, Pₙ)
            Qₙ = intLegendre(n-1; T = T, a = T(binf), b = T(bsup))
            push!(polynomials, Qₙ)
        end
        new{T}(polynomials, derivpolynomials, ordermax - ordermin + 1, ordermin, ordermax, T(binf), T(bsup))
    end
end

@inline Base.firstindex(::P1IntLegendreGenerator) = 1
@inline Base.eachindex(ilg::P1IntLegendreGenerator) = eachindex(ilg.polynomials)
@inline Base.getindex(ilg::P1IntLegendreGenerator, n::Int) =  ilg.polynomials[n] 
@inline getpolynomial(ilg::P1IntLegendreGenerator) = ilg.polynomials
@inline getderivpolynomial(ilg::P1IntLegendreGenerator) = ilg.derivpolynomials

function P1IntLegendreGenerator(mesh::Mesh, T::Type = Float64; kwargs...)
    generators = P1IntLegendreGenerator(T; kwargs...)
    size = generators.size * (lastindex(mesh) - 1)
    indices_cells       = zeros(Int, size, 1)
    indices_generators  = zeros(Int, size, 1)
    normalisation       = ones(Int, size)
    for i ∈ eachindex(mesh)[1:end-1]
        for n ∈ 1:generators.size
            indices_cells[(i-1) * generators.size + n, 1]       = i
            indices_generators[(i-1) * generators.size + n, 1]  = n
        end
    end
    PolynomialBasis(generators, mesh, size, indices_cells, indices_generators, normalisation) 
end
