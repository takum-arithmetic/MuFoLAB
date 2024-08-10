# See LICENSE file for copyright and license details.
using Suppressor
@suppress_err begin
	using MatrixDepot
end

push!(LOAD_PATH, "src/")
@suppress_err begin
	using TestMatrices
end

function is_acceptable_test_matrix(t::TestMatrix)
	if t.rank != t.m
		# reject matrices that do not have full rank
		return false
	end

	return true
end

# select those matrices from the SuiteSparse Matrix Collection with NNZ <= 50.000
predicate = sp(:) & @pred(nnz <= 50_000)

generate_test_matrices(
	"out/sparse_test_matrices.jld2",
	predicate;
	is_acceptable_test_matrix = is_acceptable_test_matrix,
)
