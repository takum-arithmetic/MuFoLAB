# See LICENSE file for copyright and license details.
using Suppressor
@suppress_err begin
	using MatrixDepot
end

push!(LOAD_PATH, "src/")
@suppress_err begin
	import TestMatricesGenerator
end

TestMatricesGenerator.generate_stochastic_test_matrices(;
	file_name = "out/stochastic_test_matrices.jld2",
	count = 1000,
	n = 100,
	density_min = 0.2,
	density_max = 0.3,
)
