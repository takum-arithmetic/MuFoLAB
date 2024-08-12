# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
import QR
import TestMatrices

function solve_qr(A::AbstractMatrix, b::AbstractVector)
	Q, R = QR.qr_givens(A)
	Q_full = Q * spdiagm(ones(typeof(A[1,1]), (size(A, 1))))
	z = Q_full' * b
	return R \ z
end

# TODO
test_matrices = TestMatrices.get_test_matrices(:sparse)
filter!(t -> (t.m == t.n && t.rank == t.m && t.nnz in 200:1000), test_matrices)

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = SolverExperimentParameters(;
				solver = solve_qr,
				preconditioner = nothing,
			),
			number_types = Experiments.all_number_types,
			test_matrices = test_matrices,
		),
	),
)
