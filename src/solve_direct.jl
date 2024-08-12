# See LICENSE file for copyright and license details.
push!(LOAD_PATH, "src/")
using Experiments
import TestMatrices

function solve_direct(A::AbstractMatrix, b::AbstractVector)
	return A \ b
end

# TODO
test_matrices = TestMatrices.get_test_matrices(:sparse)
filter!(t -> (t.m == t.n && t.rank == t.m && t.nnz in 200:1000), test_matrices)
test_matrices = test_matrices[1:1]

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = SolverExperimentParameters(;
				solver = solve_direct,
				preconditioner = nothing,
			),
			number_types = Experiments.all_number_types,
			test_matrices = test_matrices,
		),
	),
)
