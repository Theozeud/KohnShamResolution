##################################################################################
#                          Plot Piecewise LaurentPolynomial
##################################################################################
using Plots

function plot_plp(p::PiecewiseLaurentPolynomial; label = false, kwargs...) 
    plt = plot(kwargs...)
    X = []
    for i ∈ p.index
            push!(X, LinRange(p.mesh[i], p.mesh[i+1], 25)...)
    end
    if firstindex(p.mesh) ∉ p.index
        pushfirst!(X, first(p.mesh))
    end
    if lastindex(p.mesh) ∉ p.index
        push!(X, last(p.mesh))
    end
    plot!(X, p.(X), label = label)
    plt
end

function plot_plp!(plt, p::PiecewiseLaurentPolynomial; kwargs...) 
    X = []
    for i ∈ p.index
            push!(X, LinRange(p.mesh[i], p.mesh[i+1], 25)...)
    end
    if firstindex(p.mesh) ∉ p.index
        pushfirst!(X, first(p.mesh))
    end
    if lastindex(p.mesh) ∉ p.index
        push!(X, last(p.mesh))
    end
    plot!(plt, X, p.(X), kwargs...)
    plt
end