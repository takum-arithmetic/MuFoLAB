# See LICENSE file for copyright and license details.
module LU

using LinearAlgebra
using SparseArrays

function lu(A::SparseMatrixCSC{T, Int64}) where {T <: AbstractFloat}
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

	return LowerTriangular(L), UpperTriangular(U)
end

end
