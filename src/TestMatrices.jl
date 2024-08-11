# See LICENSE file for copyright and license details.
module TestMatrices

using CodecZlib
using JLD2
using LinearAlgebra
using SparseArrays
using Suppressor
@suppress_err begin
	using MatrixDepot
end

export TestMatrix, generate_test_matrices, get_test_matrices

struct TestMatrix
	name::String
	M::SparseMatrixCSC{Float64, Int64}
	m::Int
	n::Int
	nnz::Int
	rank::Int
	is_symmetric::Bool
	is_positive_definite::Bool
end

function TestMatrix(M::AbstractMatrix, name::String)
	# test if M is symmetric before type cleanup
	is_symmetric = isa(M, Symmetric)

	# clean up the matrix type
	M = SparseMatrixCSC{Float64, Int64}(M)

	# determine the rest using the clean M
	matrix_rank = rank(M)
	is_positive_definite = isposdef(M)
	matrix_nnz = nnz(M)

	return TestMatrix(
		name,
		M,
		M.m,
		M.n,
		matrix_nnz,
		matrix_rank,
		is_symmetric,
		is_positive_definite,
	)
end

function generate_test_matrices(
	file_name::String,
	predicate::Any;
	is_acceptable_test_matrix::Union{Nothing, Function} = nothing,
)
	# download the desired sparse matrices, which is necessary
	# only once, but keep it quiet
	@suppress_err begin
		MatrixDepot.load(predicate)
	end

	# generate list of names
	matrix_name_list = mdlist(predicate)

	# generate the test matrix array
	test_matrices = Vector{TestMatrix}()
	for (index, matrix_name) in enumerate(matrix_name_list)
		local t

		try
			t = TestMatrix(matrixdepot(matrix_name), matrix_name)
		catch e
			if isa(e, InexactError)
				# This exception is thrown when the input
				# matrix is, for instance, complex. However,
				# we want it to at least be an AbstractMatrix
				# of type Float64.
				continue
			else
				throw(e)
			end
		end

		# optionally use the provided filtering function
		if is_acceptable_test_matrix === nothing ||
		   is_acceptable_test_matrix(t)
			push!(test_matrices, t)
		end
	end

	# sort the test matrix array by name
	sort!(test_matrices; by = t -> t.name)

	# print the file name we are writing to standard output
	println(file_name)

	# store the test matrix array in a julia data file (JLD2)
	jldopen(file_name, "w"; compress = true) do file
		return file["test_matrices"] = test_matrices
	end
end

function get_test_matrices(type::Symbol)
	if type == :sparse
		array_file = "out/sparse_test_matrices.jld2"
	else
		throw(ArgumentError("The TestMatrix array type must be :sparse"))
	end

	return load(array_file, "test_matrices")
end

end
