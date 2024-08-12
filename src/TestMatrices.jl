# See LICENSE file for copyright and license details.
module TestMatrices

using CodecZlib
using JLD2
using LinearAlgebra
using SparseArrays
using Suppressor

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

function get_test_matrices(type::Symbol)
	if type == :sparse
		array_file = "out/sparse_test_matrices.jld2"
	else
		throw(ArgumentError("The TestMatrix array type must be :sparse"))
	end

	return load(array_file, "test_matrices")
end

end
