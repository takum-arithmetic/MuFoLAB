# See LICENSE file for copyright and license details.
module LU

using LinearAlgebra
using SparseArrays

function lu(A::SparseMatrixCSC{T, Int64}) where {T <: AbstractFloat}
	# determine full pivotisation (rows and columns) using
	# UMFPACK's lu decomposer and the matrix casted to Float64
	lud = LinearAlgebra.lu(Float64.(A))
	permutation_row = lud.p
	permutation_col = lud.q

	# we now run the LU decomposition without pivotisation on PAQ,
	# where P is the row permutation and Q is the column permutation,
	# given by construction we know that (PAQ)[j,j] != for all j in 1:n
	A = A[permutation_row, permutation_col]

	# Determine n
	n = size(A, 1)

	# Initialise L as the sparse identity matrix
	L = spdiagm(ones(T, n))

	# Initialise U as an empty sparse matrix
	U = spzeros(T, n, n)

	# Iterate over all columns
	for j in 1:n
		x =
			[L[:, 1:(j - 1)] [
				zeros(j - 1, n - j + 1)
				sparse(I, n - j + 1, n - j + 1)
			]] \ A[:, j]
		U[1:j, j] = x[1:j]
		L[(j + 1):n, j] = x[(j + 1):n] ./ x[j]
	end

	return LowerTriangular(L), UpperTriangular(U), permutation_row, permutation_col
end

end
