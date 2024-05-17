
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
exchcorr(km::KohnShamExtended) = km.exc
charge(km::KohnShamExtended) = km.z
nbelec(km::KohnShamExtended) = km.N