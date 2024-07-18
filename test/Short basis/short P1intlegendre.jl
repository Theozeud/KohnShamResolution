using KohnShamResolution
using Plots

# Parameters of the discretization
using DoubleFloats
T = Float64
Rmin = 0
Rmax = 100
Nmesh = 32
m = linmesh(Rmin,Rmax,Nmesh)
normalize = true
ordermin = 2
ordermax = 2
left = false
right = false

basis = ShortP1IntLegendreBasis(m, T; ordermin = ordermin, ordermax = ordermax, left = left, right = right, normalize = normalize, Rcut = 37)

# Plots elements of the basis
X = LinRange(-1, 1, 10000)
plt_elements = plot(legend = false)
for b ∈ basis.basisVector
    for i ∈ eachindex(b.elements)
        plot!(plt_elements, X, b.elements[i].(X))
    end
end
plt_elements

# Plots of basis
X = LinRange(Rmin, Rmax, 10000)
plt_basis = plot(legend = false)
for i ∈ 1:length(basis)
    p = build_basis(basis, i)
    plot!(plt_basis, X, p.(X))
end
plt_basis


# Plots of the derivatives of the basis
deriv_basis = deriv(basis)
X = LinRange(Rmin, Rmax, 10000)
plt_derivbasis = plot(legend = false)
for i ∈ 1:length(deriv_basis)
    p = build_basis(deriv_basis, i)
    plot!(plt_derivbasis, X, p.(X))
end
plt_derivbasis


# Assembly matrix
M₀  = mass_matrix(basis)
M₋₁ = weight_mass_matrix(basis, -1)
M₋₂ = weight_mass_matrix(basis, -2)
A   = mass_matrix(deriv(basis))

nb1 = 0
for I∈CartesianIndices(M₋₁)
    if isnan(M₋₁[I]) || abs(M₋₁[I])==Inf
        global nb1 += 1
    end
end
println("Thre is $nb1 mistake coefficiens of M₋₁.")

nb2 = 0
for I∈CartesianIndices(M₋₂)
    if isnan(M₋₂[I]) || abs(M₋₂[I])==Inf
        nb2 += 1
    end
end
println("Thre is $nb2 mistake coefficiens of M₋₂.")


function integrate2(rf::RationalFraction{T}, a::Real, b::Real; geomfun = false) where T
    @assert !haslog(rf.num) && !haslog(rf.denom)
    if iszero(rf.num) || (degmax(rf.num) == 0 && abs(rf.num[0]) < sqrt(eps(typeof(rf.num[0]))))
        NewT = promote_type(T, typeof(a), typeof(b))
        return zero(NewT)
    end
    @show simplification = iszero(diveucl2(rf.num, rf.denom)[2])
    if geomfun && !simplification
        if degmax(rf.denom) > 2
            @show "here geom 3  "
            @error "Impossibility to integrate a rational fraction with geometric functions for denominator of degree more than three."
        elseif degmax(rf.denom) == 2
            @show "here geom 2  "
            NewT = promote_type(T, typeof(a), typeof(b))
            val = zero(NewT)
            for k ∈ eachindex(rf.num)
                val += rf.num[k] * KohnShamResolution._integration_monome_over_deg2(k, rf.denom[2], rf.denom[1], rf.denom[0], a, b)
            end 
            return val
        elseif degmax(rf.denom) == 1
            @show "here geom 1  "
            NewT = promote_type(T, typeof(a), typeof(b))
            val = zero(NewT)
            for k ∈ eachindex(rf.num)
                val += rf.num[k] * KohnShamResolution._integration_monome_over_deg1(k, rf.denom[1],rf.denom[0], a, b)
            end 
            return val
        elseif degmax(rf.denom) == 0
            @show "here geom 0  "
            return integrate(rf.num, a, b)/rf.denom[0]
        end
    else
        if degmax(rf.num) ≥ degmax(rf.denom)
            @show "here 1 "
            return integrate(inEntpart2(rf), a, b)
        else
            if degmax(rf.denom) > 2
                @show "here 2 "
                return KohnShamResolution.integrate(fraction_decomp(rf), a, b)
            elseif degmax(rf.denom) == 2
                @show "here 3 "
                A = rf.num[1]
                B = rf.num[0]
                C = rf.denom[2]
                D = rf.denom[1]
                E = rf.denom[0]
                return KohnShamResolution._integrate_partial12(A, B, C, D, E, a, b) 
            elseif degmax(rf.denom) == 1
                A = rf.num[0]
                B = rf.denom[1]
                C = rf.denom[0]
                @show "here 4 "
                return KohnShamResolution._integrate_partial01(A, B, C, a, b)
            end
        end
    end
end
function inEntpart2(rf::RationalFraction)
    @assert !haslog(rf.num) && !haslog(rf.denom)
    Q,R = diveucl2(rf.num, rf.denom)
    newrf = RationalFraction(R, rf.denom)
    SommeRationalFraction(Q, [newrf])
end
function diveucl2(p::LaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where{TP, TQ}
    @assert !haslog(p) && !haslog(q)
    @assert degmin(p) ≥ 0 && degmin(q) ≥ 0
    @assert !iszero(q)
    NewT = promote_type(TP,TQ)
    if degmax(p) < degmax(q)
        return (Monomial(0, NewT(0)), p)
    end
    _q = (q/q[end])
    _q = _q /_q[end]
    Q = Monomial(0, NewT(0))
    R = KohnShamResolution.convert(NewT, p)
    while degmax(R) ≥ degmax(q)
        @show Q += Monomial(degmax(R) - degmax(q), NewT(R[end]))
        @show _q
        @show R -= _q * Monomial(degmax(R) - degmax(q), NewT(R[end]))
    end
    return (round(Q/NewT(q[end])), R)
end
Base.round(p::LaurentPolynomial; digits = 14) = LaurentPolynomial(round.(p.coeffs; digits = digits), degmin(p), haslog(p), p.coeff_log)

for I ∈ CartesianIndices(M₀)
    (ib1, iib1) = KohnShamResolution.find_basis(basis, I[1])
    (ib2, iib2) = KohnShamResolution.find_basis(basis, I[2])
    spb1 = KohnShamResolution.getbasis(basis, ib1)
    spb2 = KohnShamResolution.getbasis(basis, ib2)
    FP =  KohnShamResolution.getpolynomial(spb1, iib1, 1)
    FQ =  KohnShamResolution.getpolynomial(spb2, iib2, 1)
    Finvϕ =  KohnShamResolution.getinvshift(spb1, iib1, 1)
    Fweight_shift = Monomial(-1) ∘ Finvϕ
    Fdinvϕ = Finvϕ[1]
    spb1.elements.binf, spb1.elements.bsup
    int = integrate((FP*FQ*Fweight_shift.num)/Fweight_shift.denom, spb1.elements.binf, spb1.elements.bsup; geomfun = true)
    if isnan(int)
        println(I)
        println(FQ*FP)
        println(Fweight_shift)
        println(diveucl2(FP*FQ*Fweight_shift.num, Fweight_shift.denom))
    end
end

