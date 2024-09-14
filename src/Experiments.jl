# See LICENSE file for copyright and license details.
module Experiments

using Base.Threads
using BFloat16s
using CSV
using DataFrames
using Float128Conversions
using LinearAlgebra
using Printf
using Quadmath
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
		M = Float128.(t.M)
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


			if t.absolute_minimum < floatmin(experiment.number_types[j]) ||
			   t.absolute_maximum > floatmax(experiment.number_types[j])
				# some matrix entries are out of bounds
				measurement[j, i] = nothing
			else
				# call the main get_measurement function identified
				# by the type of parameters
				measurement[j, i] = get_measurement(
					experiment.parameters,
					experiment.number_types[j],
					M,
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
				length(experiment.test_matrices),
				length(type_names) + 1,
			),
			:auto,
		)

		# assign the first column to be the matrix names
		df[!, 1] = matrix_names

		# fill the DataFrame with values from R
		for i in 1:length(type_names)
			df[!, i + 1] = [
				if m == nothing || isnan(getfield(m, field_name))
					NaN
				else
					getfield(m, field_name)
				end for m in
				experiment_results.measurement[i, :]
			]
		end

		# Set the row names as the matrix names
		rename!(df, [Symbol("matrix\\type"); Symbol.(type_names)])

		experiment_name = chopsuffix(basename(PROGRAM_FILE), ".jl")

		# create output directory
		directory_name = "out/" * experiment_name
		if !(isdir(directory_name))
			mkdir(directory_name)
		end

		# generate file name
		file_name = directory_name * "/" * String(field_name) * ".csv"

		# write the DataFrame to the target file
		CSV.write(file_name, df)

		# print the written file name to standard output for the witness file
		println(file_name)
	end
end

# the solver experiments compare the solutions of linear systems of equations
# Ax = b via the infinity norm.
@kwdef struct SolverExperimentParameters <: AbstractExperimentParameters
	solver::Function
	preconditioner::Union{Nothing, Function}
end

struct SolverExperimentMeasurement <: AbstractExperimentMeasurement
	absolute_error::Float128
	relative_error::Float128
	logarithmic_relative_error::Float128
end

function get_logarithmic_relative_error(approx::AbstractVector, exact::AbstractVector)
	# Gustafson defines the log-rel-error as
	#
	#	relErr(approx, exact) = {
	#		NaR	     | sgn(approx) != sgn(exact)
	#		|ln(approx/exact)|   | otherwise
	#	}
	#
	# For the extension to vectors, we have to incorporate the signs
	# in a way. Given one can pull in the ||exact|| in the relative error
	# term as
	#
	#	|| (approx - exact) / ||exact|| ||,
	#
	# we also do it element-wise here: First check the signs
	# of all entries for inconsistencies; this also catches the cases
	# where either entry is zero and the other isn't.
	# If there are any, return NaN immediately.
	if sign.(approx) != sign.(exact)
		return NaN
	end

	# All entries match in sign. Compute ln(approx/exact) elementwise
	# for the other entries and return the relative error as the norm
	# of this vector. Set 0/0=1 by convention.
	zero_log(a, b) = (iszero(a) && iszero(b)) ? zero(typeof(a)) : log(a / b)

	return norm([zero_log(approx[i], exact[i]) for i in 1:length(exact)], 2)
end

function get_measurement(
	parameters::SolverExperimentParameters,
	::Type{T},
	A::SparseMatrixCSC{Float128, Int64},
) where {T <: AbstractFloat}
	local x, x_approx

	# Generate right side b such that the solution is (1,...,1).
	# We compute b in full precision and only then reduce it
	# to the target type
	x = ones(Float128, size(A, 1))
	b = A * x
	b_approx = T.(b)

	# Overwrite the matrix with its reduced form
	A_approx = T.(A)

	try
		if parameters.preconditioner != nothing
			A_approx, b_approx = parameters.preconditioner(A_approx, b_approx)
		end
		x_approx = parameters.solver(A_approx, b_approx)
	catch
		return nothing
	end

	# compute errors
	exact = x
	approx = Float128.(x_approx)
	err = exact - approx

	absolute_error = norm(err, 2)
	relative_error = norm(err, 2) / norm(exact, 2)
	logarithmic_relative_error = get_logarithmic_relative_error(approx, exact)

	return SolverExperimentMeasurement(
		absolute_error,
		relative_error,
		logarithmic_relative_error,
	)
end

end
