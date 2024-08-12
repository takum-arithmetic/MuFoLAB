# See LICENSE file for copyright and license details.
push!(LOAD_PATH, "src/")

using Suppressor
@suppress_err begin
	using MatrixDepot
	using TestMatricesGenerator
end

# select those matrices from the SuiteSparse Matrix Collection with NNZ <= 50.000
predicate = sp(:) & @pred(nnz <= 50_000)

generate_test_matrices("out/sparse_test_matrices.jld2", predicate)
