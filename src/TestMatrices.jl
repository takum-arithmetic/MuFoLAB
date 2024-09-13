# See LICENSE file for copyright and license details.
module TestMatrices

using CodecZlib
using JLD2
using LinearAlgebra
using SparseArrays
using Suppressor

struct TestMatrix
	name::String
	M::SparseMatrixCSC{Float64, Int64}
	m::Int
	n::Int
	nnz::Int
	absolute_minimum::Float64
	absolute_maximum::Float64
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

	# determine absolute matrix for absolute minimum and maximum
	M_abs = abs.(M)

	return TestMatrix(
		name,
		M,
		M.m,
		M.n,
		matrix_nnz,
		minimum(M_abs[M_abs .!= 0]),
		maximum(M_abs),
		matrix_rank,
		is_symmetric,
		is_positive_definite,
	)
end

function get_test_matrices(type::Symbol; filter_function::Union{Nothing, Function} = nothing)
	if type == :sparse
		array_file = "out/sparse_test_matrices.jld2"
	else
		throw(ArgumentError("The TestMatrix array type must be :sparse"))
	end

	# load test matrices from the file
	test_matrices = load(array_file, "test_matrices")

	# filter the test matrices
	if filter_function != nothing
		filter!(filter_function, test_matrices)
	end

	# honour the request for reduced test data
	if "--reduced-test-data" in ARGS
		# obtain matrices with reasonable size
		filter!(t -> (t.nnz in 100:2000), test_matrices)

		# get the first 200
		test_matrices = test_matrices[1:min(200, end)]
	end

	return test_matrices
end

end
