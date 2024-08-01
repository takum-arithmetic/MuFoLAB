# See LICENSE file for copyright and license details.
include("Utilities.jl")
using .Utilities

function solve_direct(A::AbstractMatrix, b::AbstractVector)
	return A \ b
end

Utilities.run_solver_experiment(; solver = solve_direct, preconditioner = nothing)
