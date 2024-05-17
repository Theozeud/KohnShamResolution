function groundstate(model::AbstractDFTModel; kwargs...)
    solver = init(model; kwargs...)
    solve!(solver; kwargs...)
    makesolution()
end


function init(model::AbstractDFTModel; kwargs...)
    # Init mesh
    # Init matrix
    # Init Cache
    # Init Initial condition
end


function solve(solver::KhonShamSolver; kwargs...)
    while #stopping criteria
        loopheader()
        performstep!()
        loopfooter()
    end
end


function makesolution()

end
