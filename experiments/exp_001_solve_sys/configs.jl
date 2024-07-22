# External Dependencies
using LinearAlgebra
using SparseArrays
using MatrixDepot
using Takums
using SoftPosit
using BFloat16s
using Plots
using BenchmarkTools
using DataFrames
using CSV


# Local Dependencies
include("../../common/infnorm.jl")
include("../../common/stopcrit.jl")


# Sparse Matrix Test Suite (see https://sparse.tamu.edu)
ld = listdata(@pred(10^1 <= n <= 2 * 10^2 && n == m && nnz / n > 10));
Ms = sort([ld[i].id for i in eachindex(ld)]);

# Number Types
Ts = [
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
];
println("") # prevent output of Ts in main script hack
