# See LICENSE file for copyright and license details.
using GraphIO.EdgeList, Graphs
using JLD2
using LinearAlgebra
using MatrixMarket
using SparseArrays

push!(LOAD_PATH, "src/")
import TestMatrices

let
	# define the relevant file names
	graph_file_list_file = "out/graph_files"
	output_file = "out/graph_test_matrices.jld2"

	# generate the test matrix array
	test_matrices = Vector{TestMatrices.TestMatrix}()

	# read the list of graph files and iterate over it
	graph_files = readlines(graph_file_list_file)
	for graph_file in graph_files
		# read in the adjacency matrix from the file depending on suffix
		A = nothing

		if endswith(graph_file, ".edges")
			# edge list file
			graph = nothing
			try
				graph = loadgraph(
					graph_file,
					"graph_key",
					EdgeListFormat(),
				)
			catch e
				println(
					stderr,
					"Failed reading $(graph_file)",
				)
				rethrow(e)
			end

			A = Float64.(adjacency_matrix(graph))
		elseif endswith(graph_file, ".mtx")
			# MatrixMarket file
			try
				A = Float64.(mmread(graph_file))
			catch e
				println(
					stderr,
					"Failed reading $(graph_file)",
				)
				rethrow(e)
			end
		else
			throw(
				ArgumentError(
					"Invalid file type in graph list",
				),
			)
		end

		# check if the adjacency matrix is non-square
		m, n = size(A)
		if m != n
			# While an adjacency matrix _must_ be symmetric and
			# anything else is amathematical, usually this is just an
			# error in the size specification and the 'excess' block
			# is empty. Check if it is.
			if m > n
				if norm(A[(n + 1):m, 1:n], 1) == 0.0
					# all good, we can crop
					A = A[1:n, 1:n]
				else
					# if we cannot crop, we just fill in
					# to obtain a mxm matrix, assuming
					# the input data intended to just
					# drop isolated nodes
					A = hcat(
						A,
						zeros(
							m,
							m -
							n,
						),
					)
				end
			else # n > m
				if norm(A[1:m, (m + 1):n], 1) == 0.0
					# all good
					A = A[1:m, 1:m]
				else
					# if we cannot crop, we just fill in
					# to obtain a mxm matrix, assuming
					# the input data intended to just
					# drop isolated nodes
					A = vcat(
						A,
						zeros(
							n -
							m,
							n,
						),
					)
				end
			end
		end

		# symmetrise the adjacency matrix using average symmetrisation
		if !issymmetric(A)
			A = (A + A') / 2.0
		end

		# compute the degree matrix, ignoring isolated nodes such
		# that they are assigned eigenvalue 1
		degree_vector = sum(A; dims = 2)[:]
		D_inv_sqrt = Diagonal(
			ifelse.(
				degree_vector .> 0,
				1.0 ./ sqrt.(abs.(degree_vector)),
				0.0,
			),
		)

		# compute the normalised Laplacian
		Lnorm = I - D_inv_sqrt * A * D_inv_sqrt

		# enforce symmetry because of possible previous rounding
		# errors
		Lnorm = triu(Lnorm, 0) + triu(Lnorm, 1)'

		# drop any zeros that might have come up
		dropzeros!(Lnorm)

		# generate a test matrix, extracting the graph ID making use
		# of the fact that the graph file name always has the
		# same format
		t = TestMatrices.TestMatrix(
			split(
				graph_file[(1 + length(
					"out/graphs/",
				)):end],
				".",
			)[1],
			Lnorm,
			size(Lnorm, 1),
			size(Lnorm, 2),
			nnz(Lnorm),
			NaN,
			NaN,
			NaN,
			0,
			issymmetric(Lnorm), # just to make sure
			false, # who cares?
		)
		push!(test_matrices, t)
	end

	# sort the test matrix array by name
	sort!(test_matrices; by = t -> t.name)

	# print the file name we are writing to standard output
	println(output_file)

	# store the test matrix array in a julia data file (JLD2)
	jldopen(output_file, "w"; compress = true) do file
		return file["test_matrices"] = test_matrices
	end
end
