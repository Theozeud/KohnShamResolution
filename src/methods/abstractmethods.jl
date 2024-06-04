abstract type AbstractKohnShamResolutionMethod end

for_sphericalsymmetry(::AbstractKohnShamResolutionMethod) = "No Information..."
for_cylindricalsymmetry(::AbstractKohnShamResolutionMethod) = "No Information..."
ismethod_for_model(::AbstractKohnShamResolutionMethod, ::AbstractDFTModel) = false

abstract type AbstractKohnShamCache end