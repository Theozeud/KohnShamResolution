function ShortP1IntLegendreBasis(mesh::OneDMesh, T::Type = Float64; normalize::Bool = false, left = false, right = false, ordermin = 2, ordermax = 2, first = true, Rcut = last(mesh))
    intleg = ShortIntLegendreBasis(mesh, T; normalize = normalize, ordermin = ordermin, ordermax = ordermax, Rcut = Rcut)
    p1 = ShortP1Basis(mesh, T; normalize = normalize, left = left, right = right)
    if first
        return CombineShortPolynomialBasis(p1, intleg)
    else
        return CombineShortPolynomialBasis(intleg, p1)
    end
end  
