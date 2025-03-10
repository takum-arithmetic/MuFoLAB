# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
using LinearAlgebra
using MicroFloatingPoints
using IterativeSolvers
using IncompleteLU
using SparseArrays
import TestMatrices

function solve_gmres_ilu(
	A::AbstractMatrix,
	b::AbstractVector,
	preparation::SolverExperimentPreparation,
)
	relative_tolerances = Dict(
		1 => eltype(A)(Float64(sqrt(eps(Floatmu{4,3})))),
		2 => eltype(A)(Float64(sqrt(eps(Float16)))),
		4 => eltype(A)(Float64(sqrt(eps(Float32)))),
		8 => eltype(A)(Float64(sqrt(eps(Float64)))),
	)

	x, history = gmres(
		A,
		b;
		Pl = ilu(A),
		reltol = relative_tolerances[sizeof(eltype(A))],
		log = true,
	)

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
					t.rank == t.m &&
					t.nnz in 1:10000
				),
			),
		),
	),
)
