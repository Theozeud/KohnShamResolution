# Modelisation of exchange correlation

abstract type ExchangeCorrelation end

struct NoExchangeCorrelation <: ExchangeCorrelation end

isthereExchangeCorrelation(::ExchangeCorrelation) = true
isthereExchangeCorrelation(::NoExchangeCorrelation) = false

# SlaterXα

struct SlaterXα <: ExchangeCorrelation end

exc(::SlaterXα, rho) = - 3/4 * (3/π)^(1/3) * rho^(4/3)
vxc(::SlaterXα, rho) = - (3/π)^(1/3) * rho^(1/3)


## LSDA MODEL

struct LSDA <: ExchangeCorrelation end

exchange_density_spin_compensated(rho)  = -3/4 * (3/π)^(1/3) * rho^(4/3)
exchange_density_spin_polarized(rho)    = -3/4 * (6/π)^(1/3) * rho^(4/3)

function exchange(::LSDA, rho_up,rho_down)
    return 1/2 * (ex_density_spin_compensated(2*rho_up) + ex_density_spin_compensated(2*rho_down))
end 

function vxc_LSDA_up(rho_up)
    return  -(6/π * rho_up )^(1/3)
end 

function vxc_LSDA_down(rho_down)
    return  -(6/π * rho_down )^(1/3)
end 


function G(rs, A, α₁, β₁, β₂, β₃, β₄, p)
    # Spin interpolation formula
    tmp = 2*A*(β₁*rs^(1/2) + β₂*rs + β₃*rs^(3/2) + β₄^(1+p))
    -2*A*(1+α₁*rs) * (log(tmp +1) - log(tmp))
end


const LSDA_CORRELATION_PARAMETERS =
      #ϵc(rs,0)   # ϵc(rs,1)   # -αc(rs)
    [     1.0           1.0        1.0;  # p 
     0.031091      0.015545   0.016887;  # A
      0.21370       0.20548    0.11125;  # α₁
       7.5957       14.1189     10.357;  # β₁
       3.5876        6.1977     3.6231;  # β₂
       1.6382        3.3662    0.88026;  # β₃
      0.49294       0.62517    0.49671   # β₄
    ]  

const LCP = LSDA_CORRELATION_PARAMETERS


function rs(rho)
    # Density parameters 
    (3/(4π*ρ))^(1/3)
end

function correlation(rho_down, rhow_up)
    rho = rho_down + rhow_up

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