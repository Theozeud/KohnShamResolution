abstract type AbstractKohnShamResolutionMethod end

for_sphericalsymmetry(::AbstractKohnShamResolutionMethod) = "No Information..."
for_cylindricalsymmetry(::AbstractKohnShamResolutionMethod) = "No Information..."


abstract type AbstractKohnShamCache end