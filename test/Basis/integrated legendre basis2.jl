using KohnShamResolution
using Plots

# Test on integrated Legendre basis
T = Float64
Rmin = 1
Rmax = 5
Nmesh = 5
m = linmesh(Rmin, Rmax, Nmesh)

ord = 3
basis = IntLegendreBasis(m, T; order = ord, left = false, right = true)

X = LinRange(-1, 1,1000)
plt = plot(legend = :outertopright)
for p ∈ basis.polynomial
    @show integrate(p*p,-1,1)
    plot!(X,p.(X))
end
plt

using UnPack
"""
    mass_matrix2(ilb::IntLegendreBasis)

Compute the mass matrix for the Integrated-Legendre Polynomials basis. 
    
The upper-right block is the mass matrix of the P1 basis. The lower-left block is a diagonal
made of the interaction of polynomials of order bigger than one, in ascending order of 
orders, arranged mesh by mesh. 

"""
function mass_matrix2(ilb::IntLegendreBasis)
    @unpack left, right, mesh, polynomial, order = ilb
    nbhf = KohnShamResolution.nb_hat_functions(ilb)
    siz  = KohnShamResolution.size(ilb)
    T = KohnShamResolution.bottom_type(ilb)
    A = zeros(T, (KohnShamResolution.size(ilb), KohnShamResolution.size(ilb)))

    #Block of hat function
    if left
        A[1,1] = (mesh[2] -  mesh[1])/3 * 3/2
    else
        A[1,1] = (mesh[3] -  mesh[1])/3 * 3/2
    end
    for I ∈ 2:nbhf-1
        A[I,I]   = (mesh[I+1 + !left] -  mesh[I-1 + !left])/3 * 3/2
        A[I,I-1] = (mesh[I+1 + !left] -  mesh[I + !left])/6   * 3/2
        A[I-1,I] = A[I,I-1]
    end
    A[nbhf, nbhf-1] = (mesh[nbhf + !left] -  mesh[nbhf - 1 + !left])/6 * 3/2
    A[nbhf-1, nbhf] = A[nbhf, nbhf-1]
    if right
        A[nbhf,KohnShamResolution.nb_hat_functions(ilb)] = (mesh[nbhf + !left] -  mesh[nbhf - 1 + !left])/3 * 3/2
    else
        A[nbhf,KohnShamResolution.nb_hat_functions(ilb)] = (mesh[nbhf + 1 + !left] -  mesh[nbhf + !left])/3 * 3/2
    end


    # Diagonal for polynomial of higer order

    idxmsh= 1
    icount = 1
    for I ∈ nbhf + 1 : siz
        A[I,I] = 2/(mesh[idxmsh + 1]-mesh[idxmsh])
        icount += 1
        if icount == order
            idxmsh += 1
            icount = 1
        end
    end

    # Interaction hat functions and polynomial of higher order
    hf_up = KohnShamResolution.getpolynomial(ilb, 1)
    hf_down = KohnShamResolution.getpolynomial(ilb, 2)

    for i ∈ eachindex(mesh)[1:end-1]
        if i == 1
            if left
                for n ∈ 2:order
                    Qₙ = KohnShamResolution.getpolynomial(ilb, n+1)
                    I = nbhf + n - 1
                    A[I, i]   = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_down, Qₙ, -1, 1)
                    A[I, i+1] = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_up, Qₙ, -1 ,1)
                end
            else
                for n ∈ 2:order
                    Qₙ = KohnShamResolution.getpolynomial(ilb, n+1)
                    I = nbhf + n - 1
                    A[I, i]   = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_up, Qₙ, -1, 1)
                end
            end
        elseif i == lastindex(mesh)-1
            for n ∈ 2:order
                Qₙ = KohnShamResolution.getpolynomial(ilb, n+1)
                I = (siz - order + 1) + n - 1
                A[I, i] = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_down, Qₙ, -1 ,1)
            end
        else 
            for n ∈ 2:order
                Qₙ = KohnShamResolution.getpolynomial(ilb, n+1)
                I = nbhf  + (order - 1) * (i-1) + n - 1
                A[I, i]   = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_down, Qₙ, -1, 1)
                A[I, i+1] = 2/(mesh[i + 1]-mesh[i]) * scalar_product(hf_up, Qₙ, -1 ,1)
            end
        end
    end
    A
end