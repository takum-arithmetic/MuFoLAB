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

function solve(
	A::AbstractMatrix,
	b::AbstractVector,
	permutation_rows::AbstractVector,
	permutation_columns::AbstractVector,
)
	# first determine the scaling matrix, whose diagonal entries are
	# the row sums of A
	C = Diagonal(1 ./ sum(A; dims = 2)[:])

	# Apply the scaling, row and column permutations to A, yielding PCAS,
	# where P is the row permutation and S is the column permutation
	# matrix
	PCAS = (C * A)[permutation_rows, permutation_columns]

	# Perform a LU decomposition without pivotisation on PCAS, which
	# works given we know that (PCAS)[j,j] != 0 holds for all j in 1:n
	# by construction
	L, U = LU.lu(PCAS)

	# The linear system is of the form Ax=b. Applying scaling and a row
	# permutation on both sides yields PCAx=PCb =: y, then multiplying
	# the identity matrix S*S^inv to the right of A yields
	# (PCAS)(S^inv*x)=y, and now we can substitute our fully
	# pivoted LU decomposition, yielding LU(S^inv*x)=y. Defining
	# z := S^inv*x (and noting that x=Sz) we obtain the system
	# LUz = y, which we solve in two stages: The outer system is
	# solved via
	#
	#	w := U*z = L\y
	#
	# and then we have left an inner system Uz = w that is solved via
	#
	#	z = U\w
	#
	y = (C * b)[permutation_rows]
	w = L \ y
	z = U \ w

	# We obtain x via x=Sz
	return z[invperm(permutation_columns)]
end

end
