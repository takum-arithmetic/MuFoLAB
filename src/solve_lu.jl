# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
using LinearAlgebra
import LU
import TestMatrices

function get_lu_row_and_column_permutations(A::SparseMatrixCSC{Float64, Int64})
	# As a preparation we determine the full pivotisation (rows and columns)
	# using UMFPACK's lu decomposer
	lud = LinearAlgebra.lu(A)

	return lud.p, lud.q
end

function solve_lu(A::AbstractMatrix, b::AbstractVector, preparation::SolverExperimentPreparation)
	# iteration count is nominally 1
	return LU.solve(A, b, preparation.permutation_rows, preparation.permutation_columns),
	1
end

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = SolverExperimentParameters(;
				get_row_and_column_permutations = get_lu_row_and_column_permutations,
				solver = solve_lu,
				preconditioner = nothing,
			),
			number_types = Experiments.all_number_types,
			test_matrices = TestMatrices.get_test_matrices(
				:sparse;
				filter_function = t -> (
					# quadratic and full rank
					t.m == t.n &&
					t.rank == t.m &&
					t.nnz in 1:10000
				),
			),
		),
	),
)
