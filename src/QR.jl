# See LICENSE file for copyright and license details.
module QR

using LinearAlgebra
using SparseArrays

export qr_givens

function qr_givens(A::SparseMatrixCSC{T, Int64}) where {T <: AbstractFloat}
	# copy A into R
	R = copy(A)

	# generate empty vector of givens rotations
	#givens_rotations = Vector{LinearAlgebra.Givens{T}}(undef, 0)
	Q = LinearAlgebra.Rotation{T}([])

	# iterate over all columns of A (using our "copy" R)
	for j in 1:(R.n)
		# iterate over all non-zero row entries from the bottom up
		for i in reverse(R[:, j].nzind)
			# we only need to process rows below the diagonal
			if i <= j
				break
			end

			# get the current column to process, as a sparse vector
			col = R[:, j]

			# all non-zero entries between (including) the diagonal
			# and (excluding) i are possible partners, those below are
			# already zero
			possible_partners =
				col.nzind[col.nzind .>= j .&& col.nzind .< i]

			if length(possible_partners) == 0
				# when we have no possible partners the
				# only possible partner is the diagonal entry
				possible_partners = [j]
			end

			# we select the best of the possible partners by
			# looking at the row R[i,j+1:end] (i.e. all entries
			# to the right of which we want to 'mirror away';
			# all entries to the left are zero by construction)
			# and finding the indices where the entries are
			# zero.
			zero_inds = setdiff(
				1:(R.n - (j + 1) + 1),
				R[i, (j + 1):(R.n)].nzind,
			)

			# determine the partner row by finding out which
			# of the partner rows has the highest number of matching
			# zero entries. By construction we get the smallest
			# index if there are multiple minima, which is
			# also good as we want the partner to be as high
			# as possible.
			#
			# NOTE: Further optimisation possible by also weighing
			# in other factors, e.g. norm of r to obtain
			# values close to 1, etc.
			partner = possible_partners[argmax(
				length.([
					intersect(
						setdiff(
							1:(R.n - (j + 1) + 1),
							R[
								p,
								(j + 1):(R.n),
							].nzind,
						),
						zero_inds,
					) for
					p in
					possible_partners
				]),
			)]

			# Generate Givens rotation of indices partner,i
			g = LinearAlgebra.givens(
				R[partner, j],
				R[i, j],
				partner,
				i,
			)
			#push!(givens_rotations, g[1])

			# Apply G to Q
			Q = Q * g[1]

			# Apply G to R
			R = g[1] * R
			R[i, j] = 0.0
			dropzeros!(R)
		end
	end

	return Q', R
end

end
