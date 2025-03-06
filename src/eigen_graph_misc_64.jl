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
				tolerance = 1e-12,
			),
			number_types = [Float64, LinearTakum64, Posit64],
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
