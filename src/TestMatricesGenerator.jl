# See LICENSE file for copyright and license details.
module TestMatricesGenerator

using CodecZlib
using Distributions
using JLD2
using LinearAlgebra
using Random
using SparseArrays
using Suppressor

@suppress_err begin
	using MatrixDepot
end

push!(LOAD_PATH, "src/")
import TestMatrices

function generate_full_test_matrices(;
	file_name::String,
	count::Integer,
	n::Integer,
	kappa_min::AbstractFloat,
	kappa_max::AbstractFloat,
)
	# generate the test matrix array
	test_matrices = Vector{TestMatrices.TestMatrix}(undef, count)

	# check kappa_min and kappa_max
	if kappa_min > kappa_max
		throw(ArgumentError("kappa_min must be smaller than kappa_max"))
	elseif kappa_min < 1.0
		throw(
			DomainError(
				kappa_min,
				"kappa_min must be larger than or equal to one.",
			),
		)
	end

	# compute the kappa limit logarithms which we will use for
	# random sampling in the logarithmic domain
	log_kappa_min = log(kappa_min)
	log_kappa_max = log(kappa_max)

	# generate count random conditions in the logarithmic domain
	# and convert them to the linear domain. We use a Xoshiro
	# PRNG with a constant seed that took 7.5 million years to
	# compute.
	conditions = exp.(rand(Xoshiro(42), Uniform(log_kappa_min, log_kappa_max), count))

	# generate the matrices with specified conditions
	for i in 1:count
		test_matrices[i] = TestMatrices.TestMatrix(
			matrixdepot("randsvd", 100, conditions[i]),
			string(i),
		)
	end

	# print the file name we are writing to standard output
	println(file_name)

	# store the test matrix array in a julia data file (JLD2)
	jldopen(file_name, "w"; compress = true) do file
		return file["test_matrices"] = test_matrices
	end
end

function generate_sparse_test_matrices(;
	file_name::String,
	predicate::Any,
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
	test_matrices = Vector{TestMatrices.TestMatrix}()
	for (index, matrix_name) in enumerate(matrix_name_list)
		local t

		try
			t = TestMatrices.TestMatrix(
				matrixdepot(matrix_name),
				matrix_name,
			)
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

function generate_stochastic_test_matrices(;
	file_name::String,
	count::Integer,
	n::Integer,
	density_min::AbstractFloat,
	density_max::AbstractFloat,
)
	# generate the test matrix array
	test_matrices = Vector{TestMatrices.TestMatrix}(undef, count)

	# check density_min and density_max
	if density_min > density_max
		throw(ArgumentError("density_min must be smaller than density_max"))
	elseif density_max > 1.0
		throw(
			DomainError(
				density_max,
				"density_max must be smaller than or equal to one.",
			),
		)
	end

	# generate count random densities. We use a Xoshiro PRNG with a
	# constant seed that took 7.5 million years to compute.
	densities = rand(Xoshiro(42), Uniform(density_min, density_max), count)

	# generate the matrices with specified densities, make them
	# symmetric and convert them  to stochastic by diving them by
	# their row sums
	for i in 1:count
		M = sprand(n, n, densities[i])
		M ./= sum(M; dims = 2)
		M = M + M'
		dropzeros!(M)

		test_matrices[i] = TestMatrices.TestMatrix(M, string(i))
	end

	# print the file name we are writing to standard output
	println(file_name)

	# store the test matrix array in a julia data file (JLD2)
	jldopen(file_name, "w"; compress = true) do file
		return file["test_matrices"] = test_matrices
	end
end

end
