# See LICENSE file for copyright and license details.
module Experiments

push!(LOAD_PATH, "src/")
using ArnoldiMethod
using Base.Threads
using BFloat16s
using Crutches
using CSV
using DataFrames
using Float128Conversions
using LinearAlgebra
using MicroFloatingPoints
import LU
import QR
using Printf
using Random
using Quadmath
using Posits
using SparseArrays
using Takums

push!(LOAD_PATH, "src/")
import TestMatrices

export AbstractExperimentParameters,
	AbstractExperimentPreparation,
	AbstractExperimentMeasurement,
	Experiment,
	ExperimentResults,
	write_experiment_results,
	SolverExperimentParameters,
	SolverExperimentPreparation,
	SolverExperimentMeasurement,
	MPIRExperimentParameters,
	MPIRExperimentPreparation,
	MPIRExperimentMeasurement,
	EigenExperimentParameters,
	EigenExperimentPreparation,
	EigenExperimentMeasurement,
	ConversionExperimentParameters,
	ConversionExperimentPreparation,
	ConversionExperimentMeasurement

all_number_types = [
	Floatmu{4, 3},
	Floatmu{5, 2},
	#Takum8,
	LinearTakum8,
	Posit8,
	#Takum16,
	LinearTakum16,
	Posit16,
	BFloat16,
	Float16,
	#Takum32,
	LinearTakum32,
	Posit32,
	Float32,
	#Takum64,
	LinearTakum64,
	Posit64,
	Float64,
]

abstract type AbstractExperimentParameters end
abstract type AbstractExperimentPreparation end
abstract type AbstractExperimentMeasurement end

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

@kwdef struct Experiment
	parameters::AbstractExperimentParameters
	number_types::Vector{DataType}
	test_matrices::Vector{TestMatrices.TestMatrix}
end

@enum MeasurementError MatrixSingular MatrixUnderOverflow

struct ExperimentResults
	experiment::Experiment
	measurement::Matrix{Union{AbstractExperimentMeasurement, MeasurementError, Missing}}
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
	measurement = Matrix{Union{AbstractExperimentMeasurement, MeasurementError, Missing}}(
		missing,
		length(experiment.number_types),
		length(experiment.test_matrices),
	)

	# this is a matrix in which the threads mark their progress
	progress = _MeasurementState.(zeros(Integer, size(measurement)))
	print_lock = ReentrantLock()

	# do all the desired measurements
	@threads for i in 1:length(experiment.test_matrices)
		t = experiment.test_matrices[i]

		# run the preparation function that computes any general
		# things that do not change with each number type
		preparation = get_preparation(experiment.parameters, t.M)

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

			if t.absolute_minimum <
			   floatmin(experiment.number_types[j]) ||
			   t.absolute_maximum >
			   floatmax(experiment.number_types[j])
				# some matrix entries are out of bounds
				measurement[j, i] =
					MatrixUnderOverflow::MeasurementError
			else
				# call the main get_measurement function identified
				# by the type of parameters, passing in the
				# prepared data
				measurement[j, i] = get_measurement(
					experiment.parameters,
					experiment.number_types[j],
					t.M,
					preparation,
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
	# the measurement matrix contains a mix of MeasurementErrors, missing and the used
	# subtype of AbstractExperimentMeasurement. The first thing we do
	# is filter out the MeasurementErrors and missings and then check if it's homogeneously
	# one type.
	local measurement_type
	local measurement_type_instance

	valid_measurements = experiment_results.measurement[typeof.(
		experiment_results.measurement
	) .<: AbstractExperimentMeasurement]

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
				if typeof(m) == MeasurementError
					if (
						m ==
						MatrixSingular::MeasurementError
					)
						-Inf
					elseif (
						m ==
						MatrixUnderOverflow::MeasurementError
					)
						Inf
					else
						throw(
							DomainError(
								m,
								"Unhandled enum type",
							),
						)
					end
				elseif isnan(
					getfield(
						m,
						field_name,
					),
				)
					# if a NaN happened otherwise it is due to
					# a singularity, we book it as an Inf
					Inf
				elseif typeof(m) == Missing
					throw(
						ErrorException(
							"There are missing measurements",
						),
					)
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
	get_row_and_column_permutations::Union{Nothing, Function}
	solver::Function
	preconditioner::Union{Nothing, Function}
end

@kwdef struct SolverExperimentPreparation <: AbstractExperimentPreparation
	x_exact::Vector{Float128}
	b_exact::Vector{Float128}
	permutation_rows::Union{Nothing, Vector{Int64}}
	permutation_columns::Union{Nothing, Vector{Int64}}
end

# according to Carson and Higham, 2017
function get_pseudorandom_x_and_b(A::SparseMatrixCSC{Float64, Int64})
	# we get the UInt64 hash of the matrix and seed a PRNG with it
	prng = Xoshiro(hash(A))

	# generate a pseudorandom right hand side b and normalise it
	# such that ||b||_inf = 1
	b_exact = rand(prng, Float128, size(A, 1))
	b_exact /= norm(b_exact, Inf)

	# determine x by solving Ax=b in quad precision using our own
	# LU implementation as Julia's backslash operator can only deal
	# with Float32 and Float64 matrices. We use UMFPACK, though, to
	# obtain our row and column permutations from the Float64 A
	qrd = LinearAlgebra.qr(A)

	x_exact = QR.solve(Float128.(A), b_exact, qrd.prow, qrd.pcol)

	return x_exact, b_exact
end

function get_preparation(parameters::SolverExperimentParameters, A::SparseMatrixCSC{Float64, Int64})
	# generate x and b with a pseudorandom b that satisfies ||b||_inf = 1
	x_exact, b_exact = get_pseudorandom_x_and_b(A)

	# get the row and column permutations via the parameter function
	local permuation_rows, permutation_columns

	if parameters.get_row_and_column_permutations != nothing
		permutation_rows, permutation_columns =
			parameters.get_row_and_column_permutations(A)
	else
		permutation_rows = nothing
		permutation_columns = nothing
	end

	return SolverExperimentPreparation(;
		x_exact = x_exact,
		b_exact = b_exact,
		permutation_rows = permutation_rows,
		permutation_columns = permutation_columns,
	)
end

struct SolverExperimentMeasurement <: AbstractExperimentMeasurement
	absolute_error::Float128
	relative_error::Float128
	logarithmic_relative_error::Float128
	iteration_count::Float128
end

function get_measurement(
	parameters::SolverExperimentParameters,
	::Type{T},
	A::SparseMatrixCSC{Float64, Int64},
	preparation::SolverExperimentPreparation,
) where {T <: AbstractFloat}
	local x_approx, iteration_count

	# reduce the exact right-hand-side b and the given matrix to
	# the target type
	b_approx = T.(preparation.b_exact)
	A_approx = T.(A)

	try
		if parameters.preconditioner != nothing
			A_approx, b_approx = parameters.preconditioner(
				A_approx,
				b_approx,
			)
		end
		x_approx, iteration_count =
			parameters.solver(A_approx, b_approx, preparation)
	catch e
		if isa(e, SingularException)
			# No problemo, we have provisioned for this
			return MatrixSingular::MeasurementError
		else
			# We just rethrow any other (unknown) exception
			# as this most likely indicates an error in the
			# code rather than a numerical phenomenon
			rethrow(e)
		end
	end

	# compute errors
	exact = preparation.x_exact
	approx = Float128.(x_approx)
	err = exact - approx

	absolute_error = norm(err, 2)
	relative_error = norm(err, 2) / norm(exact, 2)
	logarithmic_relative_error = get_logarithmic_relative_error(approx, exact)

	return SolverExperimentMeasurement(
		absolute_error,
		relative_error,
		logarithmic_relative_error,
		iteration_count,
	)
end

# the iterative refinement experiments compare the solutions of linear
# systems of equations Ax = b via the infinity norm.
@kwdef struct MPIRExperimentParameters <: AbstractExperimentParameters
	low_precision_type::DataType
	working_precision_type::DataType
	high_precision_type::DataType
	tolerance::Float64
	maximum_iteration_count::Unsigned
end

@kwdef struct MPIRExperimentPreparation <: AbstractExperimentPreparation
	x_exact::Vector{Float128}
	b_exact::Vector{Float128}
end

function get_preparation(parameters::MPIRExperimentParameters, A::SparseMatrixCSC{Float64, Int64})
	# generate x and b with a pseudorandom b that satisfies ||b||_inf = 1
	x_exact, b_exact = get_pseudorandom_x_and_b(A)

	return MPIRExperimentPreparation(; x_exact = x_exact, b_exact = b_exact)
end

struct MPIRExperimentMeasurement <: AbstractExperimentMeasurement
	absolute_error::Float128
	relative_error::Float128
	logarithmic_relative_error::Float128
	iteration_count::Float128
end

function get_measurement(
	parameters::MPIRExperimentParameters,
	::Type{T},
	A::SparseMatrixCSC{Float64, Int64},
	preparation::MPIRExperimentPreparation,
) where {T <: AbstractFloat}
	local x_approx

	# get the row and column permutations for A via UMFPACK
	lud = LinearAlgebra.lu(A)
	permutation_rows = lud.p
	permutation_columns = lud.q

	# convert the permuted exact solution b to working and high precision
	b_working = parameters.working_precision_type.(preparation.b_exact[permutation_rows])
	b_high = parameters.high_precision_type.(preparation.b_exact[permutation_rows])

	# obtain PAS, which is A with row- and column permutations
	PAS_working =
		parameters.working_precision_type.(
			A[permutation_rows, permutation_columns]
		)
	PAS_high = parameters.high_precision_type.(A[permutation_rows, permutation_columns])

	# initialise solution vector x to zero, in working precision
	x_working = zeros(parameters.working_precision_type, size(PAS_working, 2))

	# determine the residual (b - PAS * x) in high precision as it's cheap
	residual_high = b_high - PAS_high * parameters.high_precision_type.(x_working)

	# perform a LU decomposition of PAS in low precision
	L_low, U_low = LU.lu(parameters.low_precision_type.(PAS_working))

	# cast the low precision results to working precision
	L_working = parameters.working_precision_type.(L_low)
	U_working = parameters.working_precision_type.(U_low)

	iteration_count = nothing
	for i in 0:(parameters.maximum_iteration_count)
		# determine the desired upper bound on the residual
		residual_upper_bound =
			size(PAS_working, 2) *
			parameters.working_precision_type(
				parameters.tolerance,
			) *
			(
				norm(b_working, Inf) +
				norm(PAS_working, Inf) *
				norm(x_working, Inf)
			)

		# break if we get below this upper bound
		if norm(parameters.working_precision_type.(residual_high), Inf) <
		   residual_upper_bound
			iteration_count = i
			break
		end

		# otherwise, update x (and consequently the residual)
		local e_working
		try
			e_working =
				U_working \ (
					L_working \
					parameters.working_precision_type.(
						residual_high
					)
				)
		catch e
			if isa(e, SingularException)
				# No problemo, we have provisioned for this
				return MatrixSingular::MeasurementError
			else
				# We just rethrow any other (unknown) exception
				# as this most likely indicates an error in the
				# code rather than a numerical phenomenon
				rethrow(e)
			end
		end

		x_working = x_working + e_working
		residual_high =
			residual_high -
			PAS_high * parameters.high_precision_type.(e_working)
	end

	if iteration_count == nothing
		# we exceeded the maximum iteration count
		return MatrixUnderOverflow::MeasurementError
	end

	# Apply column permutation to x_working to obtain solution for A
	x_working = x_working[permutation_columns]

	# compute errors
	exact = preparation.x_exact
	approx = Float128.(x_working)
	err = exact - approx

	absolute_error = norm(err, 2)
	relative_error = norm(err, 2) / norm(exact, 2)
	logarithmic_relative_error = get_logarithmic_relative_error(approx, exact)

	return MPIRExperimentMeasurement(
		absolute_error,
		relative_error,
		logarithmic_relative_error,
		iteration_count,
	)
end

# the eigen experiments determine a set of eigenvalues for a given matrix,
# e.g. using Krylov subspace methods
@kwdef struct EigenExperimentParameters <: AbstractExperimentParameters
	which::Symbol # :LM (largest abs), :LR (most pos), :SR (most neg)
	eigenvalue_count::Int64
	tolerance::Float128
end

@kwdef struct EigenExperimentPreparation <: AbstractExperimentPreparation
	start_vector_exact::Vector{Float128}
	eigenvalues_exact::Vector{Float128}
	eigenvectors_exact::Matrix{Float128} # column vectors of matrix
end

function get_pseudorandom_v(A::SparseMatrixCSC{Float64, Int64})
	# we get the UInt64 hash of the matrix and seed a PRNG with it
	prng = Xoshiro(hash(A))

	# generate a pseudorandom vector and normalise it such that
	# ||v||_inf = 1. There is no need to normalise it, but we do
	# it anyway.
	v = rand(prng, Float128, size(A, 1))
	v /= norm(v, Inf)

	return v
end

function get_preparation(parameters::EigenExperimentParameters, A::SparseMatrixCSC{Float64, Int64})
	# check if A is symmetric
	if !issymmetric(A)
		throw(ArgumentError("Input matrix is not symmetric"))
	end

	# generate a pseudorandom Float128 start vector
	start_vector_exact = get_pseudorandom_v(A)

	# compute the exact eigenvalues and eigenvectors based on the
	# number of them requested in the parameters. The tolerance
	# is set to the limits of 128 bits.
	#
	# We use the ArnoldiMethod.jl package, as it implements the
	# Arnoldi method with Krylov-Schur restart, which is superior
	# to the ARPACK method of implictly-restarted deflation, which
	# is mathematically equivalent but more complicated to implement.
	decomposition, history = partialschur(
		Float128.(A);
		nev = parameters.eigenvalue_count,
		tol = 1e-24,
		which = parameters.which,
		v1 = start_vector_exact,
	)

	# verify that we converged and all eigenvalues are real (which
	# should be, as A is symmetric)
	if !isreal(decomposition.eigenvalues)
		throw(ArgumentError("Eigenvalues are not real"))
	end

	if !history.converged
		# We just carry on here as divergence in Float128 implies
		# divergence in any smaller types.
		println(stderr, "Eigenvalues are unconverged...")

		return EigenExperimentPreparation(;
			start_vector_exact = start_vector_exact,
			eigenvalues_exact = zeros(
				Float128,
				parameters.eigenvalue_count,
			),
			eigenvectors_exact = zeros(
				Float128,
				size(A, 1),
				parameters.eigenvalue_count,
			),
		)
	end

	# given all eigenvalues are real, we can directly obtain the
	# eigenvectors from the Schur decomposition, which are the columns
	# of Q. For convenience, we simply store the matrix and by
	# convention imply the column vectors
	#
	# We also do that because in the general case (Arnoldi instead of
	# Lanczos) the Schur decomposition yields a tridiagonal matrix
	# (not diagonal), so obtaining the eigenvectors from that is
	# a completely separate step that seems to be so scary that even
	# KrylovKit and ArnoldiMethod didn't dare to implement it themselves
	# and simply call into ARPACK. Thus we stay away from this and
	# simply look at symmetric matrices. Overall, this second step
	# is computationally cheap and the 'meat' is obtaining the
	# Schur decomposition.
	#
	# Crop them to the number of requested eigenvalues/eigenvectors,
	# as the process may yield more than we asked for.

	# The eigenvectors are the columns of Q, because A is symmetric.
	# However, eigenvectors are only unique up to their signs, so
	# we need to normalise them in some way. We do that by looking
	# at the signs of the respective first entries, and flipping the
	# respective eigenvector when the first entry is negative.
	eigenvectors_exact_firstsigns =
		(
			Float128(1.0) .-
			Float128(2.0) .* (
				decomposition.Q[
					1,
					1:(parameters.eigenvalue_count),
				] .< 0.0
			)
		)'

	return EigenExperimentPreparation(;
		start_vector_exact = start_vector_exact,
		eigenvalues_exact = decomposition.eigenvalues[1:(parameters.eigenvalue_count)],
		eigenvectors_exact = decomposition.Q[
			:,
			1:(parameters.eigenvalue_count),
		] .* eigenvectors_exact_firstsigns,
	)
end

struct EigenExperimentMeasurement <: AbstractExperimentMeasurement
	eigenvalues_absolute_error::Float128
	eigenvalues_relative_error::Float128
	eigenvalues_logarithmic_relative_error::Float128
	eigenvectors_absolute_error::Float128
	eigenvectors_relative_error::Float128
	eigenvectors_logarithmic_relative_error::Float128
	matrix_vector_product_count::Int64
end

function get_measurement(
	parameters::EigenExperimentParameters,
	::Type{T},
	A::SparseMatrixCSC{Float64, Int64},
	preparation::EigenExperimentPreparation,
) where {T <: AbstractFloat}
	local decomposition, history

	# reduce the start vector and the given matrix to the target type
	start_vector_approx = T.(preparation.start_vector_exact)
	A_approx = T.(A)

	try
		# compute the requested eigenvalues and eigenvectors
		decomposition, history = partialschur(
			A_approx;
			nev = parameters.eigenvalue_count,
			tol = T(parameters.tolerance),
			which = parameters.which,
			v1 = start_vector_approx,
		)
	catch e
		if isa(e, TaskFailedException)
			# No problemo, we have provisioned for this
			return MatrixSingular::MeasurementError
		else
			# We just rethrow any other (unknown) exception
			# as this most likely indicates an error in the
			# code rather than a numerical phenomenon
			return MatrixSingular::MeasurementError
			#rethrow(e) TODO
		end
	end

	# verify that we converged and all eigenvalues are real (which
	# should be, as A is symmetric); if we didn't converge we can
	# assume that the matrix became singular
	if !history.converged || !isreal(decomposition.eigenvalues)
		return MatrixSingular::MeasurementError
	end

	# same as in the preparation, given A is symmetric the Schur
	# decomposition yields the eigenvectors directly.
	# Crop them to the number of 'requested' eigenvalues and
	# eigenvectors, as the process may yield more.
	eigenvalues_approx = decomposition.eigenvalues[1:(parameters.eigenvalue_count)]

	# The eigenvectors are the columns of Q, because A is symmetric.
	# However, eigenvectors are only unique up to their signs, so
	# we need to normalise them in some way. We do that by looking
	# at the signs of the respective first entries, and flipping the
	# respective eigenvector when the first entry is negative.
	eigenvectors_approx_firstsigns =
		(
			T(1.0) .-
			T(2.0) .* (
				decomposition.Q[
					1,
					1:(parameters.eigenvalue_count),
				] .< 0.0
			)
		)'
	eigenvectors_approx =
		decomposition.Q[:, 1:(parameters.eigenvalue_count)] .*
		eigenvectors_approx_firstsigns

	# compute eigenvalue errors
	exact = preparation.eigenvalues_exact
	approx = Float128.(eigenvalues_approx)
	err = exact - approx

	eigenvalues_absolute_error = norm(err, 2)
	eigenvalues_relative_error = norm(err, 2) / norm(exact, 2)
	eigenvalues_logarithmic_relative_error =
		get_logarithmic_relative_error(approx, exact)

	# compute eigenvector errors
	exact = preparation.eigenvectors_exact
	approx = Float128.(eigenvectors_approx)
	err = exact - approx

	eigenvectors_absolute_error = norm(err, 2)
	eigenvectors_relative_error = norm(err, 2) / norm(exact, 2)
	eigenvectors_logarithmic_relative_error = Float128(0.0)
	#get_logarithmic_relative_error(approx, exact)

	return EigenExperimentMeasurement(
		eigenvalues_absolute_error,
		eigenvalues_relative_error,
		eigenvalues_logarithmic_relative_error,
		eigenvectors_absolute_error,
		eigenvectors_relative_error,
		eigenvectors_logarithmic_relative_error,
		history.mvproducts,
	)
end

# the conversion experiments determine the error witnessed by converting
# a given matrix to a respective type.
@kwdef struct ConversionExperimentParameters <: AbstractExperimentParameters
	# no parameters
end

@kwdef struct ConversionExperimentPreparation <: AbstractExperimentPreparation
	# no preparation
end

function get_preparation(parameters::ConversionExperimentParameters, A::SparseMatrixCSC{Float64, Int64})
	return ConversionExperimentPreparation()
end

struct ConversionExperimentMeasurement <: AbstractExperimentMeasurement
	absolute_error::Float128
	relative_error::Float128
	logarithmic_relative_error::Float128
end

function get_measurement(
	parameters::ConversionExperimentParameters,
	::Type{T},
	A::SparseMatrixCSC{Float64, Int64},
	preparation::ConversionExperimentPreparation,
) where {T <: AbstractFloat}
	# compute errors by converting the matrix to the target type,
	# expand them all as vectors for easier processing
	exact = (Float128.(A))[:]
	approx = (Float128.(T.(A)))[:]
	err = exact - approx

	absolute_error = norm(err, 2)
	relative_error = norm(err, 2) / norm(exact, 2)
	logarithmic_relative_error = get_logarithmic_relative_error(approx, exact)

	return ConversionExperimentMeasurement(
		absolute_error,
		relative_error,
		logarithmic_relative_error,
	)
end

end
