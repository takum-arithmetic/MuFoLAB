# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
using LinearAlgebra
using IterativeSolvers
using IncompleteLU
using SparseArrays
import TestMatrices

function solve_gmres_ilu(A::AbstractMatrix, b::AbstractVector, preparation::SolverExperimentPreparation)
	x, history = gmres(A, b; Pl = ilu(A), log = true)

	return x, history.iters
end

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = SolverExperimentParameters(;
				get_row_and_column_permutations = nothing,
				solver = solve_gmres_ilu,
				preconditioner = nothing,
			),
			number_types = Experiments.all_number_types,
			test_matrices = TestMatrices.get_test_matrices(
				:sparse;
				filter_function = t -> (
					# quadratic and full rank
					t.m == t.n &&
					t.rank == t.m
				),
			),
		),
	),
)
