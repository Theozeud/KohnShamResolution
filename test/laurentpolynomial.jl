using KohnShamResolution
using Test
using Plots


# Test "Piecewise LaurentPolynomial"
T = Float16
m = mesh(1:5)

hf1 = HatFunctionP1(m, 1, Float16)
hf2 = HatFunctionP1(m, 2, Float16)
hf3 = HatFunctionP1(m, 3, Float16)

@test hf3(m[2]) == 0
@test hf3(m[3]) == 1
@test hf3(m[4]) == 0
@test hf3(m[3]) isa T

lpb = LaurentPolynomialBasis([hf1,hf2,hf3])
mass_matrix(lpb, m[begin], m[end])

hfbasis = HatBasis(m, Float16)
@time mass_matrix(hfbasis, m[begin], m[end])

deriv!(hfbasis)

mass_matrix(hfbasis, m[begin], m[end])