# See LICENSE file for copyright and license details.
using BFloat16s
using LinearAlgebra
using MatrixDepot
using SoftPosit
using SparseArrays
using Takums

include("Utilities.jl")
using .Utilities

"""
    Experiment 1: Compare solutions for different types of numbers for 
    solving linear systems of equations, Ax = b. In this experiment, 
    we compare the infinity norm of the difference between the true 
    solution and the computed solution for different types of numbers. 
    The types of numbers we compare are Takums, SoftPosits, BFloat16s, and 
    Floats.  We use the SuiteSparse Matrix Collection to test the different 
    number types.  The integer ids are the 'official' ident numbers 
    assigned by the collection. Currently id ∈ 1:3000.
"""
struct SolveData
	absolute_error::Float64
end

function solve_direct(M::AbstractMatrix, t::Type{T}) where {T <: AbstractFloat}
	local x

	A = Utilities.sparsecast(t, M)
	y = ones(t, size(A, 1))
	b = A * y

	try
		x = A \ b
	catch
		x = fill(NaN, size(A, 1))
	end

	absolute_error = norm(Float64.(y) - Float64.(x), Inf)

	return SolveData(absolute_error)
end

types = [
	Takum8,
	Posit8,
	Takum16,
	Posit16_1,
	Posit16,
	BFloat16,
	Float16,
	Takum32,
	Posit32,
	Float32,
	Takum64,
	Float64,
]
matrix_indices = sort(
	getfield.(listdata(@pred(10^1 <= n <= 2 * 10^2 && n == m && nnz / n > 10)), :id),
)

results = Utilities.run_experiment(solve_direct, types, matrix_indices)

Utilities.write_experiment_csv(
	"out/solve_direct.csv",
	types,
	matrix_indices,
	getfield.(results, :absolute_error),
)
