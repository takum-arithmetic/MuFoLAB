# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
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
				tolerance = 1e-8,
			),
			number_types = [Float32, Takum32, Posit32],
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
