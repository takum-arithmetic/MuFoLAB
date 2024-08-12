# See LICENSE file for copyright and license details.
module Utilities

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

number_types = [
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

struct SolveMeasurement
	absolute_error::Float64
end

function run_experiment(
	f::Function,
	types::Vector{DataType},
	test_matrices::Vector{TestMatrices.TestMatrix},
	parameters::Any = nothing,
)
	# get return type of passed function and generate output matrix
	f_return_types = Base.return_types(f)
	if (length(f_return_types) != 1)
		throw("Experiment function has inconsistent return types")
	end
	results = Array{f_return_types[1], 2}(undef, length(types), length(test_matrices))

	@threads for i in 1:length(test_matrices)
		t = test_matrices[i]
		@threads for j in 1:length(types)
			println(
				stderr,
				"Started $(t.name) for type $(String(nameof(types[j])))",
			)
			results[j, i] = f(
				types[j].(t.M),
				types[j],
				parameters,
			)
			println(
				stderr,
				"Finished $(t.name) for type $(String(nameof(types[j])))",
			)
		end
	end

	return results
end

function write_experiment_csv(
	file_name::String,
	types::Vector{DataType},
	test_matrices::Vector{TestMatrices.TestMatrix},
	R::AbstractMatrix,
)
	# generate array of strings
	type_names = String.(nameof.(types))
	matrix_names = (t -> t.name).(test_matrices)

	# generate CSV
	df = DataFrame(
		Matrix{typeof(R[1, 1])}(
			undef,
			length(type_names),
			length(test_matrices) + 1,
		),
		:auto,
	)

	# assign the first column to be the type names
	df[!, 1] = type_names

	# fill the DataFrame with values from R
	for i in 1:length(test_matrices)
		df[!, i + 1] = R[:, i]
	end

	# Set the column names as the matrix names
	rename!(df, [Symbol("type\\matrix"); Symbol.(matrix_names)])

	# write the DataFrame to the target file
	CSV.write(file_name, df)

	# print the written file name to standard output for the witness file
	return println(file_name)
end

# the solver experiments compare the solutions of linear systems of equations
# Ax = b via the infinity norm.
struct SolverParameters
	solver::Function
	preconditioner::Union{Nothing, Function}
end

function solver_experiment_callback(
	A::SparseMatrixCSC{T, Int64},
	t::Type{T},
	parameters::SolverParameters,
) where {T <: AbstractFloat}
	local x

	y = ones(t, size(A, 1))
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

	return SolveMeasurement(absolute_error)
end

function run_solver_experiment(; solver::Function, preconditioner::Union{Function, Nothing})
	test_matrices = TestMatrices.get_test_matrices(:sparse)
	filter!(t -> (t.m == t.n && t.rank == t.m && t.nnz in 200:1000), test_matrices)
	test_matrices = test_matrices[1:1]

	results = run_experiment(
		solver_experiment_callback,
		Utilities.number_types,
		test_matrices,
		SolverParameters(solver, preconditioner),
	)

	experiment_name = chopsuffix(basename(PROGRAM_FILE), ".jl")

	return write_experiment_csv(
		"out/" * experiment_name * "-absolute_error.csv",
		number_types,
		test_matrices,
		getfield.(results, :absolute_error),
	)
end

end
