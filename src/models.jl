# Modelisation of exchange correlation

abstract type ExchangeCorrelation end

struct NoExchangeCorrelation <: ExchangeCorrelation end

isthereExchangeCorrelation(::ExchangeCorrelation) = true
isthereExchangeCorrelation(::NoExchangeCorrelation) = false

# SlaterXα MODEL

struct SlaterXα <: ExchangeCorrelation end

exc(::SlaterXα, rho) = - 3/4 * (3/π)^(1/3) * rho^(4/3)
vxc(::SlaterXα, rho) = - (3/π)^(1/3) * rho^(1/3)


## LSDA MODEL

struct LSDA <: ExchangeCorrelation end

exc(lsda::LSDA, rhoUP, rhoDOWN) = ex(lsda, rhoUP,rhoDOWN) + ec(lsda, rhoUP,rhoDOWN)
vxcUP(lsda::LSDA, rhoUP, rhoDOWN) = vxUP(lsda, rhoUP) + vcUP(lsda, rhoUP, rhoDOWN)
vxcDOWN(lsda::LSDA, rhoUP, rhoDOWN) = vxDOWN(lsda, rhoDOWN) + vcDOWN(lsda, rhoUP, rhoDOWN)

ex(::LSDA, rhoUP,rhoDOWN) =  -3/4 * (6/π)^(1/3) * (rhoUP^(4/3) + rhoDOWN^(4/3))
vxUP(::LSDA, rhoUP) = -(6/π * rhoUP )^(1/3)
vxDOWN(::LSDA, rhoDOWN) = -(6/π * rhoDOWN )^(1/3)


function G(p, rs, A, α₁, β₁, β₂, β₃, β₄)
    # Spin interpolation formula
    tmp = 2*A*(β₁*rs^(1/2) + β₂*rs + β₃*rs^(3/2) + β₄^(1+p))
    -2*A*(1+α₁*rs) * (log(tmp +1) - log(tmp))
end

function ∂G(p, rs, A, α₁, β₁, β₂, β₃, β₄)
    rs12    = rs^(1/2)
    rs32    = rs * rs^(1/2)
    rsp     = rs^p
    Q0      = -2*A*(1+α₁*rs)
    Q1      = 2*A*(β₁*rs12 + β₂*rs + β₃*rs32 + β₄*rsp * rs)
    derivQ1 = A*(β₁*rs^(-1/2) + 2*β₂ + 3*β₃*rs12 + 2*(p+1)*β₄ *rsp)
    return -2 * A * α₁ * (log(Q1+1) -log(Q1)) - (Q0*derivQ1)/(Q1^2 + Q1)
end


const LSDA_CORRELATION_PARAMETERS =
      # εcPW0       # εcPW1      # -αc
    [     1.0           1.0        1.0;  # p 
     0.031091      0.015545   0.016887;  # A
      0.21370       0.20548    0.11125;  # α₁
       7.5957       14.1189     10.357;  # β₁
       3.5876        6.1977     3.6231;  # β₂
       1.6382        3.3662    0.88026;  # β₃
      0.49294       0.62517    0.49671   # β₄
    ]  

f(ξ) = ((1 + ξ)^(4/3) + (1 - ξ)^(4/3) - 2)/(2^(4/3) - 2)
derivf(ξ) = 4/3 * ((1 + ξ)^(1/3) + (1 - ξ)^(1/3))/(2^(4/3) - 2)

const invd2f0  = (9 * (2^(1/3) -1))/4                       # 1/f''(0)

function ec(::LSDA, rhoUP, rhoDOWN)
    rho = rhoDOWN + rhoUP                                   # Density
    ξ   = (rhoUP - rhoDOWN)/rho                             # Relative spin polarization
    rs  = (3/(4π*rho))^(1/3)                                # Density parameter
    return rho * εcPW(rs, ξ)
end

function εcPW(rs, ξ)
    εcPW0 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)     # Correlation energy densioty for ξ = 0
    εcPW1 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,2]...)     # Correlation energy densioty for ξ = 1
    αc    = -G(rs, LSDA_CORRELATION_PARAMETERS[:,3]...)     # Spin stiffness
    fξ    = f(ξ)                                            # f(ξ)
    ξ4    = ξ^4                                             # ξ^4
    return εcPW0 + αc * (1 - ξ4) * fξ * invd2f0 + (εcPW1 - εcPW0) * ξ4 * fξ
end

function vcUP(::LSDA, rhoUP, rhoDOWN)
    rho = rhoDOWN + rhoUP                                   # Density
    ξ   = (rhoUP - rhoDOWN)/rho                             # Relative spin polarization
    rs  = (3/(4π*rho))^(1/3)                                # Density parameter
    εcPW0 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)     # Correlation energy densioty for ξ = 0
    εcPW1 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,2]...)     # Correlation energy densioty for ξ = 1
    αc    = -G(rs, LSDA_CORRELATION_PARAMETERS[:,3]...)     # Spin stiffness
    εcPWξ  =  εcPW(rs, ξ)
    ∂rs∂εcPW0 =  ∂G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)
    ∂rs∂εcPW1 =  ∂G(rs, LSDA_CORRELATION_PARAMETERS[:,2]...)
    dαc       = -∂G(rs, LSDA_CORRELATION_PARAMETERS[:,3]...)
    fξ    = f(ξ)                                            # f(ξ)
    ξ4    = ξ^4                                             # ξ4
    ∂rs∂εcPW =  ∂rs∂εcPW0 + (∂rs∂εcPW1 - ∂rs∂εcPW0) * fξ * ξ4 + dαc * f(ξ) * (1 - ξ4) * invd2f0
    ΔεcPW10  = εcPW1 - εcPW0
    ∂ξ∂εcPW  = 4ξ^3 *fξ * (ΔεcPW10 - αc * invd2f0) + derivf(ξ) * (ξ4 * ΔεcPW10 + (1 - ξ4) * αc * invd2f0)
    return εcPWξ - rs/3 * ∂rs∂εcPW - (ξ - 1) * ∂ξ∂εcPW
end

function vcDOWN(::LSDA, rhoUP, rhoDOWN)
    rho = rhoDOWN + rhoUP                                   # Density
    ξ   = (rhoUP - rhoDOWN)/rho                             # Relative spin polarization
    rs  = (3/(4π*rho))^(1/3)                                # Density parameter
    εcPW0 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)     # Correlation energy densioty for ξ = 0
    εcPW1 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,2]...)     # Correlation energy densioty for ξ = 1
    αc    = -G(rs, LSDA_CORRELATION_PARAMETERS[:,3]...)     # Spin stiffness
    εcPWξ  =  εcPW(rs, ξ)
    ∂rs∂εcPW0 =  ∂G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)
    ∂rs∂εcPW1 =  ∂G(rs, LSDA_CORRELATION_PARAMETERS[:,2]...)
    dαc       = -∂G(rs, LSDA_CORRELATION_PARAMETERS[:,3]...)
    fξ    = f(ξ)                                            # f(ξ)
    ξ4    = ξ^4 
    ∂rs∂εcPW = ∂rs∂εcPW0 * (1 - fξ * ξ4) + ∂rs∂εcPW1 * fξ * ξ4 + dαc * f(ξ) * (1 - ξ^4) * invd2f0
    ΔεcPW10  = εcPW1 - εcPW0
    ∂ξ∂εcPW  = 4ξ^3 *fξ * (ΔεcPW10 - αc * invd2f0) + derivf(ξ) * (ξ4 * ΔεcPW10 + (1 - ξ4) * αc * invd2f0)
    return εcPWξ - rs/3 * ∂rs∂εcPW - (ξ + 1) * ∂ξ∂εcPW
end

# KohnSham Model

abstract type AbstractDFTModel end 

struct KohnShamExtended{TEXCH <: ExchangeCorrelation, TZ <: Real, TN <: Real} <: AbstractDFTModel
    z::TZ
    N::TN
    exc::TEXCH

    function KohnShamExtended(;z, N, exc::ExchangeCorrelation = NoExchangeCorrelation())
        new{typeof(exc), eltype(z), eltype(N)}(z,N,exc)
    end
end

isthereExchangeCorrelation(km::KohnShamExtended) = isthereExchangeCorrelation(km.exc)

ReducedHartreeFock(z::Real, N::Int) = KohnShamExtended(z = z, N = N)
SlaterXα(z::Real, N::Int) = KohnShamExtended(z = z, N = N, exc = SlaterXα())
LDSA(z::Real, N::Int) = KohnShamExtended(z = z, N = N, exc = LDSA())