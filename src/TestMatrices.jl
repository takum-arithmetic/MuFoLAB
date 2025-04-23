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

function graph_get_type_from_name(name::String)
	category_map = Dict(
		:graph_biological => [
			"bio/",     # Biological Networks
			"eco/",     # Ecology
			"protein/", # Proteins
			"bn/",      # Brain Networks
		],
		:graph_social => [
			"ca/",	   # Collaboration Networks
			"cit/",	   # Citation Networks
			"dynamic/",	   # Interaction/Recommendation Networks
			"econ/",	   # Economic Networks
			"email/",	   # E-Mail Networks
			"ia/",	   # Interaction Networks
			"proximity/",	   # Human Contact Network/Interaction
			"rec/",	   # Recommender Networks
			"retweet_graphs/", # Retweet Networks
			"rt/",	   # Retweet Networks
			"soc/",	   # Social Networks
			"socfb/",	   # Social Networks (Facebook)
			"tscc/",	   # Temporal Reachability Graphs
		],
		:graph_infrastructure => [
			"inf/",     # Infrastructure Networks
			"massive/", # Large-Scale Graphs (Internet)
			"power/",   # Power Grids
			"road/",    # Road Networks
			"tech/",    # Technological Networks
			"web/",     # Web Graphs
		],
		:graph_misc => [
			"dimacs/",   # DIMACS Challenge Mix
			"dimacs10/", # DIMACS10 Challenge Mix
			"graph500/", # GRAPH500 Mix
			"heter/",    # Heterogeneous Graphs
			"labeled/",  # Labeled Graphs Mix
			"misc/",     # Miscellaneous Graphs
			"rand/",     # Random Graphs
			"sc/",       # Scientific Computing
		],
	)

	for (graph_type, categories) in category_map
		for category in categories
			if startswith(name, category)
				return graph_type
			end
		end
	end

	throw(ArgumentError("Graph name $(name) could not be assigned a category"))
end

function get_test_matrices(type::Symbol; filter_function::Union{Nothing, Function} = nothing)
	if type == :full
		array_file = "out/full_test_matrices.jld2"
	elseif type == :graph_biological ||
	       type == :graph_social ||
	       type == :graph_infrastructure ||
	       type == :graph_misc
		array_file = "out/graph_test_matrices.jld2"
	elseif type == :sparse
		array_file = "out/sparse_test_matrices.jld2"
	elseif type == :stochastic
		array_file = "out/stochastic_test_matrices.jld2"
	else
		throw(
			ArgumentError(
				"The TestMatrix array type must be :full, :sparse or :stochastic",
			),
		)
	end

	# load test matrices from the file
	test_matrices = load(array_file, "test_matrices")

	# filter out the graphs that we want using our classification function
	if type in [ :graph_biological :graph_infrastructure :graph_social :graph_misc ]
		filter!(
			tm -> (graph_get_type_from_name(tm.name) == type),
			test_matrices,
		)
	end

	# filter the test matrices
	if filter_function != nothing
		filter!(filter_function, test_matrices)
	end

	# honour the request for reduced test data
	if "--reduced-test-data" in ARGS
		if type in [ :full :graph_biological :graph_infrastructure :graph_social :graph_misc ]
			# obtain matrices with reasonable size
			filter!(t -> (t.nnz in 1:10000), test_matrices)

			# get the first 50
			test_matrices = test_matrices[1:min(50, end)]
		elseif type == :sparse || type == :stochastic
			# obtain matrices with reasonable size
			filter!(t -> (t.nnz in 1:10000), test_matrices)

			# get the first 50
			test_matrices = test_matrices[1:min(50, end)]
		else
			throw(ArgumentError("Unknown matrix category $(type)"))
		end
	end

	return test_matrices
end

end
