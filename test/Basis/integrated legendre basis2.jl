using KohnShamResolution
using Plots

# Test on integrated Legendre basis
T = Float64
Rmin = 1
Rmax = 5
Nmesh = 5
m = linmesh(Rmin, Rmax, Nmesh)

ord = 1
basis = IntLegendreBasis(m, T; order = ord, left = true, right = true)

X = LinRange(-1, 1,1000)
plt = plot(legend = :outertopright)
for p ∈ basis.polynomial
    @show integrate(p*p,-1,1)
    plot!(X,p.(X))
end
plt

using UnPack
function mass_matrix2(ilb::IntLegendreBasis)
    @unpack left, right, mesh, polynomial, order = ilb

    T = KohnShamResolution.bottom_type(ilb)
    A = zeros(T, (KohnShamResolution.size(ilb), KohnShamResolution.size(ilb)))

    #Block of hat function
    if left
        A[1,1] = (mesh[2] -  mesh[1])/3
    else
        A[1,1] = (mesh[3] -  mesh[1])/3
    end
    for I ∈ 2:KohnShamResolution.nb_hat_functions(ilb)-1
        A[I,I]   = (mesh[I+1 + !left] -  mesh[I-1 + !left])/3
        A[I,I-1] = (mesh[I+1 + !left] -  mesh[I + !left])/6
        A[I-1,I] = A[I,I-1]
    end
    A[KohnShamResolution.nb_hat_functions(ilb),KohnShamResolution.nb_hat_functions(ilb)-1] = (mesh[KohnShamResolution.nb_hat_functions(ilb) + !left] -  mesh[KohnShamResolution.nb_hat_functions(ilb) - 1 + !left])/6
    A[KohnShamResolution.nb_hat_functions(ilb)-1,KohnShamResolution.nb_hat_functions(ilb)] = A[KohnShamResolution.nb_hat_functions(ilb),KohnShamResolution.nb_hat_functions(ilb)-1]
    if right
        A[KohnShamResolution.nb_hat_functions(ilb),KohnShamResolution.nb_hat_functions(ilb)] = (mesh[KohnShamResolution.nb_hat_functions(ilb) + !left] -  mesh[KohnShamResolution.nb_hat_functions(ilb) - 1 + !left])/3
    else
        A[KohnShamResolution.nb_hat_functions(ilb),KohnShamResolution.nb_hat_functions(ilb)] = (mesh[KohnShamResolution.nb_hat_functions(ilb) + 1 + !left] -  mesh[KohnShamResolution.nb_hat_functions(ilb) + !left])/3
    end


    # Diagonal for polynomial of higer order
    idxmsh= 1
    icount = 1
    for I ∈ KohnShamResolution.nb_hat_functions(ilb)+1:KohnShamResolution.size(ilb)
        A[I,I] = 2/(mesh[idxmsh + 1]-mesh[idxmsh])
        icount += 1
        if icount == order
            idxmsh += 1
            icount = 1
        end
    end

    # Interaction hat functions and polynomial of higher order
    for i ∈ eachindex(mesh)[1:end-1]
        for n ∈ 2:order
            Qₙ = KohnShamResolution.getpolynomial(ilb, n+1)
            hf_up = KohnShamResolution.getpolynomial(ilb, 1)
            hf_down = KohnShamResolution.getpolynomial(ilb, 2)
            A[KohnShamResolution.nb_hat_functions(ilb)+n-1, i]   = 2/(mesh[i + 1]-mesh[i]) * integrate(hf_up, Qₙ, -1, 1)
            A[KohnShamResolution.nb_hat_functions(ilb)+n-1, i+1] = 2/(mesh[i + 1]-mesh[i]) * integrate(hf_down, Qₙ, -1 ,1)
        end
    end
    A
end