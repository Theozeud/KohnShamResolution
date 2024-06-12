function eigen_value_OneElectronAtom(z::Real, N::Int)
    [-z/(2*n^2) for n ∈ 1:N]
end

function eigen_value_Hydrogen(N::Int)
    eigen_value_OneElectronAtom(1,N)
end

