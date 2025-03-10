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
			parameters = ConversionExperimentParameters(),
			number_types = Experiments.all_number_types,
			test_matrices = TestMatrices.get_test_matrices(
				:sparse;
				filter_function = t -> (
					t.nnz <= 50000
				),
			),
		),
	),
)
