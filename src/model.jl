abstract type AbstractDFTModel end

```
    KohnSham Model
```
struct KohnShamExtended{TZ,TN,TEXCH} <: AbstractDFTModel
    z::TZ
    N::TN
    exc::TEXCH
end

Base.size(km::KohnShamExtended) = (km.z,km.N)
@inline exchcorr(km::KohnShamExtended) = km.exc
@inline charge(km::KohnShamExtended) = km.z
@inline nbelec(km::KohnShamExtended) = km.N