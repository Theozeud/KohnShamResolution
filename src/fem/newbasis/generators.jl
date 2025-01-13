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
        polynomials = Vector{LaurentPolynomial{T}}(undef,ordermax-ordermin+1)
        derivpolynomials = Vector{LaurentPolynomial{T}}(undef,ordermax-ordermin+1) 
        for n ∈ ordermin:ordermax
            Pₙ = Legendre(n-1; T = T, a = T(binf), b = T(bsup))
            derivpolynomials[n-ordermin+1] = Pₙ
            Qₙ = intLegendre(n-1; T = T, a = T(binf), b = T(bsup))
            polynomials[n-ordermin+1] = Qₙ
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
    cells_to_indices = zeros(Int, length(mesh)-1, generators.size) 
    for i ∈ eachindex(mesh)[1:end-1]
        cells_to_indices[i,:] = ((i-1) * generators.size + 1):(i * generators.size)
        for n ∈ 1:generators.size
            indices_cells[(i-1) * generators.size + n, 1]       = i
            indices_generators[(i-1) * generators.size + n, 1]  = n
        end
    end
    PolynomialBasis(generators, mesh, size, indices_cells, indices_generators, cells_to_indices, 
                    normalisation) 
end


# P1 + Integrated Legendre Generators

struct P1IntLegendreGenerator{T} <: AbstractGenerator{T}
    polynomials::Vector{LaurentPolynomial{T}}
    derivpolynomials::Vector{LaurentPolynomial{T}}
    size::Int
    ordermax::Int
    binf::T
    bsup::T

    function P1IntLegendreGenerator(T::Type = Float64; ordermax = 2, binf::Real = -T(1), bsup::Real = T(1))
        @assert ordermax ≥ 1
        polynomials = Vector{LaurentPolynomial{T}}(undef, ordermax+1)
        derivpolynomials = Vector{LaurentPolynomial{T}}(undef, ordermax+1)
        polynomials[1] = Polynomial([one(T),one(T)], 0)
        polynomials[2] = Polynomial([one(T),-one(T)], 0)
        derivpolynomials[1] = Polynomial([one(T)], 0)
        derivpolynomials[2] = Polynomial([-one(T)], 0)
        for n ∈ 2:ordermax
            Pₙ = Legendre(n-1; T = T, a = T(binf), b = T(bsup))
            derivpolynomials[n+1] = Pₙ
            Qₙ = intLegendre(n-1; T = T, a = T(binf), b = T(bsup))
            polynomials[n+1] = Qₙ
        end
        new{T}(polynomials, derivpolynomials, ordermax + 1, ordermax, T(binf), T(bsup))
    end
end

@inline Base.firstindex(::P1IntLegendreGenerator) = 1
@inline Base.eachindex(p1ilg::P1IntLegendreGenerator) = eachindex(p1ilg.polynomials)
@inline Base.getindex(p1ilg::P1IntLegendreGenerator, n::Int) =  p1ilg.polynomials[n] 
@inline getpolynomial(p1ilg::P1IntLegendreGenerator) = p1ilg.polynomials
@inline getderivpolynomial(p1ilg::P1IntLegendreGenerator) = p1ilg.derivpolynomials

function P1IntLegendreGenerator(mesh::Mesh, T::Type = Float64; kwargs...)
    generators = P1IntLegendreGenerator(T; kwargs...)
    size = (generators.ordermax - 1)* (lastindex(mesh) - 1) +  (lastindex(mesh) - 2)
    indices_cells       = zeros(Int, size, 2)
    indices_generators  = zeros(Int, size, 2)
    normalisation       = ones(Int, size)
    cells_to_indices = zeros(Int, length(mesh)-1, generators.size)
    for i ∈ eachindex(mesh)[1:end-2]
        indices_cells[i,:]        = [i,i+1] 
        indices_generators[i,:]   = [1,2]
    end
    for i ∈ eachindex(mesh)[1:end-1]
        if i != firstindex(mesh) && i != lastindex(mesh)-1
            cells_to_indices[i,:] = ((i-1) * generators.size):(i * generators.size-1)
        elseif i == firstindex(mesh)
            cells_to_indices[1,1:end-1] = 1:(generators.size-1)
        else
            cells_to_indices[end,1:end-1] = (size - generators.size + 2):size
        end
        for n ∈ 2:generators.ordermax
            index = (lastindex(mesh) - 2) + (i-1) * (generators.size-2) + n -1
            indices_cells[index, 1]       = i
            indices_generators[index, 1]  = n+1 
        end
    end
    cells_to_indices[1,2] = 0
    cells_to_indices[end,1] = 0
    PolynomialBasis(generators, mesh, size, indices_cells, indices_generators, cells_to_indices, 
        normalisation) 
end
