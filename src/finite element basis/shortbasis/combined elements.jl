function ShortP1IntLegendreBasis(mesh::Mesh, T::Type = Float64; left = false, right = false, ordermin = 2, ordermax = 2, first = true, Rcut = last(mesh))
    if ordermax == 1
        return ShortP1Basis(mesh, T; left = left, right = right)
    else
        intleg = ShortIntLegendreBasis(mesh, T; ordermin = ordermin, ordermax = ordermax, Rcut = Rcut)
        p1 = ShortP1Basis(mesh, T; left = left, right = right)
        if first
            return CombineShortPolynomialBasis(p1, intleg)
        else
            return CombineShortPolynomialBasis(intleg, p1)
        end
    end
end  
