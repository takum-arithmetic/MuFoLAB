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
			number_types = [ Float64 ],
			test_matrices = TestMatrices.get_test_matrices(
				:graph_social;
				filter_function = t -> (
					t.n >= 12
				),
			),
		),
	),
)
