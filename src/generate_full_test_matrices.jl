# See LICENSE file for copyright and license details.
using Suppressor

push!(LOAD_PATH, "src/")
@suppress_err begin
	import TestMatricesGenerator
end

TestMatricesGenerator.generate_full_test_matrices(;
	file_name = "out/full_test_matrices.jld2",
	count = 1000,
	n = 100,
	kappa_min = 1.0,
	kappa_max = 1e15,
)
