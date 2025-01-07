
struct DFTProblem{TMo <: AbstractDFTModel, TD <: KohnShamDiscretization, TMe <: SCFMethod}
    model::TMo
    discretization::TD
    method::TMe
end

@inline model(p::DFTProblem)            = p.model
@inline discretization(p::DFTProblem)   = p.discretization
@inline method(p::DFTProblem)           = p.method
