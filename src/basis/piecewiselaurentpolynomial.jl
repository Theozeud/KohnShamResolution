struct PiecewiseLaurentPolynomial{T,TM}
    laurentpolynomials::Vector{LaurentPolynomial{T}}
    intervals::Vector{Tuple{TM,TM}}
end

