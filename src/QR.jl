# See LICENSE file for copyright and license details.
module QR

using LinearAlgebra
using Metis
using SparseArrays

function qr_givens(A::SparseMatrixCSC{T, Int64}) where {T <: AbstractFloat}
	# determine fill-reducing permutation (rows and columns) using
	# Metis and the trivial symmetrisation A+A'
	permutation, inverse_permutation = Metis.permutation(A + A')

	# copy permuted A into R
	R = copy(A[permutation, permutation])

	# generate empty vector of givens rotations
	Q = LinearAlgebra.Rotation{T}([])

	# iterate over all columns of the permuted A (using our "copy" R)
	for j in 1:(R.n)
		# non-zero entries in the column below the diagonal
		nonzero_indices = filter(ind -> (ind > j), R[:, j].nzind)

		# work from the bottom up, cancelling the nonzero entry
		# using a Givens rotation with the next nonzero entry
		# above it, or using the diagonal for the topmost entry
		for m in reverse(1:length(nonzero_indices))
			top_index    = (m == 1) ? j : nonzero_indices[m - 1]
			bottom_index = nonzero_indices[m]

			# Generate Givens rotation
			g = LinearAlgebra.givens(
				R[top_index, j],
				R[bottom_index, j],
				top_index,
				bottom_index,
			)

			# Apply G to Q
			Q = Q * g[1]

			# Apply G to R
			R = g[1] * R
			R[bottom_index, j] = 0.0
		end
	end

	dropzeros!(R)

	return Q', UpperTriangular(R), permutation, inverse_permutation
end

end
