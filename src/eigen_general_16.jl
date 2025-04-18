# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using BFloat16s
using Crutches
using Experiments
using LinearAlgebra
using Posits
using Takums
import TestMatrices

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = EigenExperimentParameters(;
				which = :LR,
				eigenvalue_count = 10,
				eigenvalue_buffer_count = 2,
				tolerance = 1e-4,
			),
			number_types = [
				BFloat16,
				Float16,
				LinearTakum16,
				Posit16,
			],
			test_matrices = TestMatrices.get_test_matrices(
				:sparse;
				filter_function = t -> (
					# quadratic and symmetric, n >= 15
					t.m == t.n &&
					t.n >= 12 &&
					t.nnz <= 20000 &&
					t.is_symmetric
				),
			),
		),
	),
)
