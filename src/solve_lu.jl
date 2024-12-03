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
	# first determine the scaling matrix, whose diagonal entries are
	# the row sums of A
	C = Diagonal(1 ./ sum(A; dims=2)[:])

	# Apply the scaling, row and column permutations to A, yielding PCAS,
	# where P is the row permutation and S is the column permutation
	# matrix
	PCAS = (C * A)[preparation.permutation_rows, preparation.permutation_columns]

	# Perform a LU decomposition without pivotisation on PCAS, which
	# works given we know that (PCAS)[j,j] != 0 holds for all j in 1:n
	# by construction
	L, U = LU.lu(PCAS)

	# The linear system is of the form Ax=b. Applying scaling and a row
	# permutation on both sides yields PCAx=PCb =: y, then multiplying
	# the identity matrix S*S^inv to the right of A yields
	# (PCAS)(S^inv*x)=y, and now we can substitute our fully
	# pivoted LU decomposition, yielding LU(S^inv*x)=y. Defining
	# z := S^inv*x (and noting that x=Sz) we obtain the system
	# LUz = y, which we solve in two stages: The outer system is
	# solved via
	#
	#	w := U*z = L\y
	#
	# and then we have left an inner system Uz = w that is solved via
	#
	#	z = U\w
	#
	y = (C * b)[preparation.permutation_rows]
	w = L \ y
	z = U \ w

	# We obtain x via x=Sz
	return z[preparation.permutation_columns], 1
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
					t.rank == t.m
				),
			),
		),
	),
)
