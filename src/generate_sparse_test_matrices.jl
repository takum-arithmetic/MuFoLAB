# See LICENSE file for copyright and license details.
using Suppressor
@suppress_err begin
	using MatrixDepot
end

push!(LOAD_PATH, "src/")
@suppress_err begin
	import TestMatricesGenerator
end

# select those matrices from the SuiteSparse Matrix Collection with NNZ <= 50.000
predicate = sp(:) & @pred(nnz <= 50_000)

TestMatricesGenerator.generate_sparse_test_matrices(;
	file_name = "out/sparse_test_matrices.jld2",
	predicate = predicate,
)
