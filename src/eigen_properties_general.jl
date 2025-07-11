# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
using LinearAlgebra
import TestMatrices

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = PropertiesExperimentParameters(),
			number_types = [Float64],
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
