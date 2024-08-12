module Experiments

using Base.Threads
using BFloat16s
using CSV
using DataFrames
using LinearAlgebra
using SoftPosit
using SparseArrays
using Takums

push!(LOAD_PATH, "src/")
import TestMatrices

export AbstractExperimentParameters,
	AbstractExperimentMeasurement,
	Experiment,
	ExperimentResults,
	write_experiment_results,
	SolverExperimentParameters,
	SolverExperimentMeasurement

all_number_types = [
	Takum8,
	Posit8,
	Takum16,
	Posit16_1,
	Posit16,
	BFloat16,
	Float16,
	Takum32,
	Posit32,
	Float32,
	Takum64,
	Float64,
]

abstract type AbstractExperimentParameters end
abstract type AbstractExperimentMeasurement end

@kwdef struct Experiment
	parameters::AbstractExperimentParameters
	number_types::Vector{DataType}
	test_matrices::Vector{TestMatrices.TestMatrix}
end

struct ExperimentResults
	experiment::Experiment
	measurement::Matrix{AbstractExperimentMeasurement}
end

function ExperimentResults(experiment::Experiment)
	results = Matrix{AbstractExperimentMeasurement}(
		undef,
		length(experiment.number_types),
		length(experiment.test_matrices),
	)

	@threads for i in 1:length(experiment.test_matrices)
		t = experiment.test_matrices[i]
		@threads for j in 1:length(experiment.number_types)
			println(
				stderr,
				"Started $(t.name) for type $(String(nameof(experiment.number_types[j])))",
			)

			# call the main get_measurement function identified
			# by the type of parameters
			results[j, i] = get_measurement(
				experiment.parameters,
				experiment.number_types[j].(t.M),
			)
			println(
				stderr,
				"Finished $(t.name) for type $(String(nameof(experiment.number_types[j])))",
			)
		end
	end

	return ExperimentResults(experiment, results)
end

function write_experiment_results(experiment_results::ExperimentResults)
	for field_name in fieldnames(typeof(experiment_results.measurement[1, 1]))
		# get underlying experiment
		experiment = experiment_results.experiment

		# generate array of strings
		type_names = String.(nameof.(experiment.number_types))
		matrix_names = (t -> t.name).(experiment.test_matrices)

		# generate CSV
		df = DataFrame(
			Matrix{
				typeof(
					getfield(
						experiment_results.measurement[
							1,
							1,
						],
						field_name,
					),
				),
			}(
				undef,
				length(type_names),
				length(experiment.test_matrices) + 1,
			),
			:auto,
		)

		# assign the first column to be the type names
		df[!, 1] = type_names

		# fill the DataFrame with values from R
		for i in 1:length(experiment.test_matrices)
			df[!, i + 1] =
				getfield.(
					experiment_results.measurement[
						:,
						i,
					],
					field_name,
				)
		end

		# Set the column names as the matrix names
		rename!(df, [Symbol("type\\matrix"); Symbol.(matrix_names)])

		# generate file name
		experiment_name = chopsuffix(basename(PROGRAM_FILE), ".jl")
		file_name =
			"out/" *
			experiment_name *
			"-" *
			String(field_name) *
			".csv"

		# write the DataFrame to the target file
		CSV.write(file_name, df)

		# print the written file name to standard output for the witness file
		return println(file_name)
	end
end

# the solver experiments compare the solutions of linear systems of equations
# Ax = b via the infinity norm.
@kwdef struct SolverExperimentParameters <: AbstractExperimentParameters
	solver::Function
	preconditioner::Union{Nothing, Function}
end

struct SolverExperimentMeasurement <: AbstractExperimentMeasurement
	absolute_error::Float64
end

function get_measurement(
	parameters::SolverExperimentParameters,
	A::SparseMatrixCSC{T, Int64},
) where {T <: AbstractFloat}
	local x

	y = ones(T, size(A, 1))
	b = A * y

	try
		if parameters.preconditioner != nothing
			A, b = parameters.preconditioner(A, b)
		end
		x = parameters.solver(A, b)
	catch
		x = fill(NaN, size(A, 1))
	end

	absolute_error = norm(Float64.(y) - Float64.(x), Inf)

	return SolverExperimentMeasurement(absolute_error)
end

end
