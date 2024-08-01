# See LICENSE file for copyright and license details.
module Utilities

using Base.Threads
using CSV
using DataFrames
using LinearAlgebra
using SparseArrays
using Suppressor

@suppress_err begin
	using MatrixDepot
end

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

function sparsecast(t::Type{T}, M::AbstractMatrix) where {T <: AbstractFloat}
	if typeof(M) <: SparseMatrixCSC
		return SparseMatrixCSC{t, Int64}(M)
	elseif typeof(M) <: Symmetric
		return Symmetric{t, SparseMatrixCSC{t, Int64}}(M)
	else
		throw("Unhandled sparse matrix type")
	end
end

function run_experiment(f::Function, types::Vector{DataType}, matrix_indices::Vector{Int64})
	# get return type of passed function and generate output matrix
	f_return_types = Base.return_types(f)
	if (length(f_return_types) != 1)
		throw("Experiment function has inconsistent return types")
	end
	results = Array{f_return_types[1], 2}(undef, length(types), length(matrix_indices))

	@threads for i in 1:length(matrix_indices)
		M = matrixdepot(sp(matrix_indices[i]))
		@threads for j in 1:length(types)
			results[j, i] = f(M, types[j])
		end
	end

	return results
end

function write_experiment_csv(
	file_name::String,
	types::Vector{DataType},
	matrix_indices::Vector{Int64},
	R::AbstractMatrix,
)
	# generate array of strings
	type_names = String.(nameof.(types))
	matrix_names = mdlist(sp(matrix_indices))

	# generate CSV
	df = DataFrame(
		Matrix{typeof(R[1, 1])}(
			undef,
			length(type_names),
			length(matrix_indices) + 1,
		),
		:auto,
	)

	# assign the first column to be the type names
	df[!, 1] = type_names

	# fill the DataFrame with values from R
	for i in 1:length(matrix_indices)
		df[!, i + 1] = R[:, i]
	end

	# Set the column names as the matrix names
	rename!(df, [Symbol("type\\matrix"); Symbol.(matrix_names)])

	# write the DataFrame to the target file
	CSV.write(file_name, df)

	# print the written file name to standard output for the witness file
	return println(file_name)
end
end
