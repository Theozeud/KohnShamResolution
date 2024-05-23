abstract type AbstractDFTModel end

```
    KohnSham Model
```
struct KohnShamExtended{TZ,TN,TEXCH,TPOT} <: AbstractDFTModel
    z::TZ
    N::TN
    exc::TEXCH
    potential::TPOT
    function KohnShamExtended(z::Real,N::Int, exc::Base.Callable = NothingFunction(), potential::Base.Callable = NothingFunction())
        new{typeof(z), typeof(N), typeof(exc), typeof(potential)}(z,N,exc,potential)
    end
end

Base.size(km::KohnShamExtended) = (km.z,km.N)
@inline exchcorr(km::KohnShamExtended) = km.exc
@inline charge(km::KohnShamExtended) = km.z
@inline nbelec(km::KohnShamExtended) = km.N
@inline potential(km::KohnShamExtended) = km.potential