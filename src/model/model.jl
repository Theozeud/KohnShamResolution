abstract type AbstractDFTModel end

```
    KohnSham Model
```
struct KohnShamExtended{TZ,TN,TEXCH,TPOT} <: AbstractDFTModel
    z::TZ
    N::TN
    exc::TEXCH = NothingFunction()
    potential::TPOT = NothingFunction()
end

Base.size(km::KohnShamExtended) = (km.z,km.N)
@inline exchcorr(km::KohnShamExtended) = km.exc
@inline charge(km::KohnShamExtended) = km.z
@inline nbelec(km::KohnShamExtended) = km.N
@inline potential(km::KohnShamExtended) = km.potential