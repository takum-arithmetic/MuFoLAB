# See LICENSE file for copyright and license details.
module QR

using LinearAlgebra
using SparseArrays

export HouseholderReflection, qr_householder, qr_givens

@kwdef struct HouseholderReflection
	normal_vector::AbstractVector{<:AbstractFloat}
	indices::Vector{Int}
end

# We overload the constructor with this function, that takes a full
# column and the index of the diagonal element to yield a householder
# reflection that maps all the non-zero entries below the diagonal to zero
function get_householder_reflection_and_lambda(;
	column::AbstractVector{T},
	diagonal_index::Int,
) where {T <: AbstractFloat}
	# identify the non-zero indices in the column and remove all which
	# are smaller than the diagonal index (we won't touch them)
	indices = filter!(i -> (i > diagonal_index), findall(!iszero, column))

	# re-add the diagonal index
	prepend!(indices, diagonal_index)

	# catch the special case where indices is empty (all other cases
	# work seamlessly). Also, just for convenience, catch the
	# case where length(indices) == 1 and indices[1] == diagonal_index,
	# because then the reflection would be trivial, too.
	if isempty(indices) || (length(indices) == 1 && indices[1] == diagonal_index)
		return nothing, zero(T)
	end

	# filter out the target vector from the column
	target = column[indices]

	# determine lambda, which is simply the norm of the target vector
	# multiplied with the sign of -target[1] to avoid numerical
	# cancellation (by convention we say sign(0)=1)
	lambda = (-target[1] >= 0) ? norm(target, 2) : -norm(target, 2)

	# The normal vector is the normalised result of
	# (target - (lambda,0,...,0)).
	raw_normal_vector = target
	raw_normal_vector[1] -= lambda
	normal_vector = raw_normal_vector ./ norm(raw_normal_vector, 2)

	# Return the reflection and lambda
	return HouseholderReflection(; normal_vector = normal_vector, indices = indices),
	lambda
end

function Base.:*(
	reflection::HouseholderReflection,
	A::AbstractVecOrMat{T},
) where {T <: AbstractFloat}
	M = copy(A)
	M[reflection.indices, :] -=
		(T(2.0) .* reflection.normal_vector) *
		(reflection.normal_vector' * M[reflection.indices, :])
	return M
end

function Base.:*(Q::Vector{HouseholderReflection}, A::AbstractVecOrMat{T}) where {T <: AbstractFloat}
	#
	# It holds
	#
	#    Q = reflection_1 * ... * reflection_n
	#
	# and thus
	#
	#    Q * A = reflection_1 * ... * reflection_n * A
	#
	# so we first apply reflection n from the left, then n-1, etc.
	#
	M = copy(A)
	for i in reverse(1:length(Q))
		M = Q[i] * M
	end
	return M
end

# given we only work with real numbers this is fine
Base.adjoint(reflection::HouseholderReflection) = reflection

#
# It holds
#
#    Q = reflection_1 * ... * reflection_n
#
# and, making use of the fact that Householder matrices are symmetric,
#
#    Q^T = (reflection_1 * ... * reflection_n)^T
#        = reflection_n^T * ... * reflection_1^T
#        = reflection_n * ... * reflection_1
#
Base.adjoint(Q::Vector{HouseholderReflection}) = reverse(Q)

function qr_householder(A::AbstractMatrix{T}) where {T <: AbstractFloat}
	# copy A into R
	R = copy(A)

	# generate empty vector of Householder reflections
	Q = HouseholderReflection[]

	# iterate over all columns of A (using our copy R)
	for column_index in 1:size(R, 2)
		# generate a householder reflection from the column vector
		# and the diagonal index, which is equal to the index of
		# the current column
		reflection, lambda = get_householder_reflection_and_lambda(;
			column = R[:, column_index],
			diagonal_index = column_index,
		)

		if reflection == nothing
			# no need for a reflection in this column, carry on
			continue
		else
			# set the column's diagonal entry to lambda and
			# the rest of the column to zero, making use of
			# the indices provided by the householder
			# transformation
			R[reflection.indices, column_index] .= zero(T)
			R[column_index, column_index] = lambda

			# apply the reflection to R (given all entries to
			# the left are zero we only have to apply it to
			# the columns to the right
			R[:, (column_index + 1):end] =
				reflection *
				R[:, (column_index + 1):end]

			# append the reflection to Q
			push!(Q, reflection)   ### ADD A SINGLE VECTOR
		end
	end

	return Q, R
end

function qr_givens(A::SparseMatrixCSC{T, Int64}) where {T <: AbstractFloat}
	# copy A into R
	R = copy(A)

	# generate empty vector of Givens rotations
	Q = LinearAlgebra.Rotation{T}([])

	# iterate over all columns of A (using our copy R)
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

	return Q', UpperTriangular(R)
end

end
