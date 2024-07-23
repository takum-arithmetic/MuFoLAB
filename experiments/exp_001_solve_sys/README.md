# Experiment 1: 

Compare the solutions solving linear systems of equations, $Ax = b$, for various floating-point number systems.
In all comparisons, the actual solution is ${\bf x} = (1,1,\dots,1)^\top$.  The infinity norm between the exact 
and computed solutions is reported.

## Run Experiment

Navigate to 

```
├── experiments/
│   ├── exp_001_solve_sys/
│   │   ├── xgesol.jl
│   │   ├── configs.jl
│   │   ├── solvesys.jl
│   │   └── README.md
```

1. Run the configurations, `configs.jl`.
2. Main experiement is in `solvesys.jl`.  Note: `@timev precomp = f.([905,906])` precompiles main function.


### Numeric Types Tested

```julia
# Numeric Types
[
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
 	Float64
] 
```

### Test-Suite

We use [SuiteSparse Matrix Collection](https://sparse.tamu.edu) (formerly the University of Florida Matrix Collection).
The collection currently assigns the `id` where `id ∈ 1:3000`.  The matrices are sparse.

#### Numeric Identifiers

Access using the numerical identifier `id` as:

```julia
using MatrixDepot
A = Array(T.(matrixdepot(sp(i))))  # where T <: AbstractFloat
```

To work with a sparse matrix in type `T`, cast as:

```
A = SparseMatrixCSC{T, Int64}(A) 
```

#### Filter Matrices of Interest

One way to filter matrices is to use `@pred`.  For example, square matrices of size 100 to 2000 with the number 
of nonzero (nnz) entries greater than a threshold (e.g., $10n$) is:

```
listdata(@pred(10^2 <= n <= 2*10^3 && n == m && nnz / n > 10 ))
```



### External Dependencies

```julia
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
```


## References 

* Weijian Zhang and Nicholas J. Higham,
  "Matrix Depot: An Extensible Test Matrix Collection for Julia," *PeerJ Comput. Sci.*, 2:e58 (2016),
  [pdf](https://peerj.com/articles/cs-58/)
