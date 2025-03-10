# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using BFloat16s
using Crutches
using Experiments
using LinearAlgebra
using MicroFloatingPoints
using Posits
using Takums
import TestMatrices

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = EigenExperimentParameters(;
				which = :LR,
				eigenvalue_count = 10,
				tolerance = 1e-2,
			),
			number_types = [Floatmu{4,3}, Floatmu{5,2}, LinearTakum8, Posit8],
			test_matrices = TestMatrices.get_test_matrices(
				:graph_misc;
				filter_function = t -> (
					t.nnz >= 100 &&
					t.nnz <= 10000
				),
			),
		),
	),
)
