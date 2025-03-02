

mutable struct CombinedMethod{typeMethod1 <: SCFMethod, typeMethod2 <: SCFMethod} <: SCFMethod
    method1::typeMethod1
    method2::typeMethod2
    state::Int
    iter_switch::Int
end


function performstep!(solver::KohnShamSolver, method::CombinedMethod)

    if solver.iter == method.iter_switch
        method.state = 2
    end

    if method.state == 1
        performstep!(solver, method.method1)
    else
        performstep!(solver, method.method2)
    end
end