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
	condition::Float64
	rank::Int
	is_symmetric::Bool
	is_positive_definite::Bool
end

function TestMatrix(M::AbstractMatrix, name::String)
	# clean up the matrix type
	if typeof(M) != Matrix{Float64}
		M = SparseMatrixCSC{Float64, Int64}(M)
	end

	m = size(M, 1)
	n = size(M, 2)

	# determine the rank
	matrix_rank = rank(M)

	# determine the condition number for square full-rank matrices
	matrix_condition = Inf
	if m == n && matrix_rank == n
		if typeof(M) != Matrix{Float64}
			# we need to pivotise manually as the SparseMatrix
			# implementation is garbage
			lud = lu(M)
			matrix_condition = cond(M[lud.p, lud.q], 1)
		else
			matrix_condition = cond(M, 1)
		end
	end

	# determine the rest using the clean M
	is_symmetric = issymmetric(M)
	is_positive_definite = isposdef(M)
	if typeof(M) != Matrix{Float64}
		matrix_nnz = nnz(M)
	else
		matrix_nnz = m * n
	end

	# determine absolute matrix for absolute minimum and maximum
	M_abs = abs.(M)

	return TestMatrix(
		name,
		M,
		m,
		n,
		matrix_nnz,
		minimum(M_abs[M_abs .!= 0]),
		maximum(M_abs),
		matrix_condition,
		matrix_rank,
		is_symmetric,
		is_positive_definite,
	)
end

function get_test_matrices(type::Symbol; filter_function::Union{Nothing, Function} = nothing)
	if type == :full
		array_file = "out/full_test_matrices.jld2"
	elseif type == :sparse
		array_file = "out/sparse_test_matrices.jld2"
	else
		throw(ArgumentError("The TestMatrix array type must be :full or :sparse"))
	end

	# load test matrices from the file
	test_matrices = load(array_file, "test_matrices")

	# filter the test matrices
	if filter_function != nothing
		filter!(filter_function, test_matrices)
	end

	# honour the request for reduced test data
	if "--reduced-test-data" in ARGS
		if type == :full
			# get the first 200
			test_matrices = test_matrices[1:min(200, end)]
		elseif type == :sparse
			# obtain matrices with reasonable size
			filter!(t -> (t.nnz in 100:2000), test_matrices)

			# get the first 200
			test_matrices = test_matrices[1:min(200, end)]
		end
	end

	return test_matrices
end

end
