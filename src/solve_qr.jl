# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
import QR
import TestMatrices

function solve_qr(A::AbstractMatrix, b::AbstractVector)
	Q, R, permutation, inverse_permutation = QR.qr(A)
	Q_full = Q * spdiagm(ones(typeof(A[1, 1]), (size(A, 1))))
	z = Q_full' * b[permutation]
	return (R \ z)[inverse_permutation]
end

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = SolverExperimentParameters(;
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
