```
    Exchange Correlation Model
```
abstract type AbstractExchangeCorrelation end

struct NoExchangeCorrelation <: AbstractExchangeCorrelation end

isthereExchangeCorrelation(::AbstractExchangeCorrelation) = true
isthereExchangeCorrelation(::NoExchangeCorrelation) = false

struct ExchangeCorrelation <: AbstractExchangeCorrelation
    exc
    vxc
end

build_SlaterXα() = ExchangeCorrelation(exc_SlaterXα, vxc_SlaterXα)
exc_SlaterXα(ρ) = -3/4 * (3/π)^(1/3) * ρ^(4/3)
vxc_SlaterXα(ρ) = - (3/π)^(1/3) * ρ^(1/3)


abstract type AbstractDFTModel end

```
    KohnSham Model
```
Base.@kwdef struct KohnShamExtended{TZ,TN,TEXCH <: ExchangeCorrelation,TPOT} <: AbstractDFTModel
    z::TZ
    N::TN
    exc::TEXCH = NoExchangeCorrelation()
    potential::TPOT = NothingFunction()
end

Base.size(km::KohnShamExtended) = (km.z,km.N)
@inline exchcorr(km::KohnShamExtended) = km.exc
@inline charge(km::KohnShamExtended) = km.z
@inline nbelec(km::KohnShamExtended) = km.N
@inline potential(km::KohnShamExtended) = km.potential


ReducedHartreeFock(z::Real, N::Int, potential::Base.Callable = NothingFunction()) = KohnShamExtended(z = z, N = N, potential = potential)
SlaterXα(z::Real, N::Int, potential::Base.Callable = NothingFunction()) = KohnShamExtended(z = z, N = N, exc = build_SlaterXα(), potential = potential)

