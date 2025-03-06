# See LICENSE file for copyright and license details.
using CSV
using DataFrames

push!(LOAD_PATH, "src/")
import TestMatrices

# get condition numbers, filter and sort them
test_matrices = TestMatrices.get_test_matrices(:sparse);
filter!(t -> (t.m == t.n && t.rank == t.m), test_matrices);
filter!(t -> (t.nnz in 1:20000), test_matrices);
conditions = sort((t -> t.condition).(test_matrices))
perc = (collect(0:(length(conditions) - 1))) ./ (length(conditions) - 1)

# write output file
df = DataFrame(; percent = perc, condition_number = conditions)
CSV.write("out/sparse_matrix_condition_numbers.csv", df)

# print output witness
println("out/sparse_matrix_condition_numbers.csv")
