# See LICENSE file for copyright and license details.
module Experiments

using Base.Threads
using BFloat16s
using CSV
using DataFrames
using LinearAlgebra
using Printf
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
	measurement::Matrix{Union{Nothing, AbstractExperimentMeasurement}}
end

@enum _MeasurementState pending = 0 processing done

let
	# export the function name globally
	global _print_progress

	# we are within the let..end scope, this simulates a static variable
	last_output_length = 0

	function _print_progress(
		progress::Matrix{_MeasurementState},
		number_types::Vector{DataType},
		test_matrices::Vector{TestMatrices.TestMatrix},
	)
		# make a local copy
		progress = copy(progress)

		# determine the number of completed and total measurements
		num_done = length(progress[progress .== done])
		num_total = length(progress)

		# generate a matrix of measurement identifiers ("matrix_name[type_name]")
		identifiers = [
			m.name * "[" * String(nameof(t)) * "]" for t in number_types,
			m in test_matrices
		]

		# get a list of identifiers that are currently active
		currently_processing_identifiers =
			identifiers[progress .== processing]

		# generate output string
		output = @sprintf "(%03i/%03i)" num_done num_total
		if length(currently_processing_identifiers) > 0
			output *= " currently processing:"
			for id in currently_processing_identifiers
				output *= " " * id
			end
		end

		# move back carriage return, write sufficiently many spaces to
		# cover the previous output, then another carriage return,
		# then the output string
		print(stderr, "\r" * (" "^last_output_length) * "\r" * output)

		# set the last output length for the next iteration
		return last_output_length = length(output)
	end
end

function ExperimentResults(experiment::Experiment)
	# this is where we store the measurements
	measurement = Matrix{Union{Nothing, AbstractExperimentMeasurement}}(
		nothing,
		length(experiment.number_types),
		length(experiment.test_matrices),
	)

	# this is a matrix in which the threads mark their progress
	progress = _MeasurementState.(zeros(Integer, size(measurement)))
	print_lock = ReentrantLock()

	# do all the desired measurements
	@threads for i in 1:length(experiment.test_matrices)
		t = experiment.test_matrices[i]
		@threads for j in 1:length(experiment.number_types)
			# set the current problem to active
			progress[j, i] = processing

			# print the current status when the lock is not
			# held, otherwise just keep going
			if trylock(print_lock)
				_print_progress(
					progress,
					experiment.number_types,
					experiment.test_matrices,
				)
				unlock(print_lock)
			end

			# cast the test matrix to the test type and check
			# if any entries are Inf or NaN
			test_matrix = experiment.number_types[j].(t.M)

			if !all(isfinite, test_matrix)
				measurement[j, i] = nothing
			else
				# call the main get_measurement function identified
				# by the type of parameters
				measurement[j, i] = get_measurement(
					experiment.parameters,
					experiment.number_types[j].(
						t.M
					),
				)
			end

			# set the current problem to done
			progress[j, i] = done
			if trylock(print_lock)
				_print_progress(
					progress,
					experiment.number_types,
					experiment.test_matrices,
				)
				unlock(print_lock)
			end
		end
	end

	# print a new line
	print(stderr, "\n")

	return ExperimentResults(experiment, measurement)
end

function write_experiment_results(experiment_results::ExperimentResults)
	# the measurement matrix contains a mix of "nothing" and the used
	# subtype of AbstractExperimentMeasurement. The first thing we do
	# is filter out the nothing and then check if it's homogeneously
	# one type.
	local measurement_type
	local measurement_type_instance

	valid_measurements =
		experiment_results.measurement[experiment_results.measurement .!= nothing]

	if length(valid_measurements) == 0
		# we only have invalid measurements
		throw(ArgumentError("All measurements are invalid"))
	else
		# check if all valid measurements are of one type
		measurement_types = union(typeof.(valid_measurements))

		if length(measurement_types) != 1
			throw(
				ArgumentError(
					"Measurement types are not homogeneous",
				),
			)
		else
			measurement_type = measurement_types[1]
			measurement_type_instance = valid_measurements[1]
		end
	end

	for field_name in fieldnames(measurement_type)
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
						measurement_type_instance,
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
			df[!, i + 1] = [
				if m == nothing
					NaN
				else
					getfield(m, field_name)
				end for m in
				experiment_results.measurement[:, i]
			]
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
		return nothing
	end

	absolute_error = norm(Float64.(y) - Float64.(x), Inf)

	return SolverExperimentMeasurement(absolute_error)
end

end
