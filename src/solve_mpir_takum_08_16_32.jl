# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
using Crutches
using Experiments
using LinearAlgebra
using Takums
import LU
import TestMatrices

write_experiment_results(
	ExperimentResults(
		Experiment(;
			parameters = MPIRExperimentParameters(;
				low_precision_type = Takum8,
				working_precision_type = Takum16,
				high_precision_type = Takum32,
				tolerance = 1.0e-03,
				maximum_iteration_count = 100,
			),
			number_types = [Float64],
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
