using KohnShamResolution
using LinearAlgebra


T = Rational{Int}

# Tests integrales on elements

intLeg = IntLegendreElements(T; ordermin = 2, ordermax = 10, binf = -T(1), bsup = T(1))

for n ∈ eachindex(intLeg)
    println("Leg deg $(n)")
    @show getderivpolynomial(intLeg,n)
    println("IntLeg deg $(n+1)")
    @show intLeg[n]
end

element_massmatrix = zeros(T, intLeg.size, intLeg.size)
for I ∈ CartesianIndices(element_massmatrix)
    element_massmatrix[I] = scalar_product(intLeg[I[1]], intLeg[I[2]], intLeg.binf, intLeg.bsup)
end

element_stiffnessmatrix = zeros(T, intLeg.size, intLeg.size)
for I ∈ CartesianIndices(element_stiffnessmatrix)
    element_stiffnessmatrix[I] = scalar_product(getderivpolynomial(intLeg,I[1]), getderivpolynomial(intLeg,I[2]), intLeg.binf, intLeg.bsup)
end

# Study 