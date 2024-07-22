"""
    Returns the maximum absolute value of each column of a matrix.
"""
function maxcols(A::Matrix{T}) where {T <: AbstractFloat}
	return y = vec(maximum(abs, A; dims = 1))::Vector{T}
end
