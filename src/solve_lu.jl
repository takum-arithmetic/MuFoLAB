# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
import LU
import TestMatrices

function solve_lu(A::AbstractMatrix, b::AbstractVector)
	L, U, permutation_row, permutation_col = LU.lu(A)

	# the linear system is of the form Ax=b. Applying a row
	# permutation on both sides yields PAx=Pb =: y, then multiplying
	# the identity matrix Q*Q^inv to the right of A yields
	# (PAQ)(Q^inv*x)=y, and now we can substitute our fully
	# pivoted LU decomposition, yielding LU(Q^inv*x). Defining
	# z := Q^inv*x (and noting that x=Qz) we obtain the system
	# LUz = y, which we solve in two stages: The outer system is
	# solved via
	#
	#	w := U*z = L\y
	#
	# and then we have left an inner system Uz = w that is solved via
	#
	#	z = U\w
	#
	y = b[permutation_row]
	w = L \ y
	z = U \ w

	# We obtain x via x=Qz
	return z[permutation_col]
end

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = SolverExperimentParameters(;
				solver = solve_lu,
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
