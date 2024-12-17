# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
using LinearAlgebra
using SparseArrays
import QR
import TestMatrices

function get_qr_row_and_column_permutations(A::SparseMatrixCSC{Float64, Int64})
	# as a preparation we determine the fill-reducing row
	# and column permutations using SPQR
	qrd = LinearAlgebra.qr(A)

	return qrd.prow, qrd.pcol
end

function solve_qr(A::AbstractMatrix, b::AbstractVector, preparation::SolverExperimentPreparation)
	# iteration count is nominally 1
	return QR.solve(A, b, preparation.permutation_rows, preparation.permutation_columns),
	1
end

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = SolverExperimentParameters(;
				get_row_and_column_permutations = get_qr_row_and_column_permutations,
				solver = solve_qr,
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
