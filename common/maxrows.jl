"""
    Returns the maximum absolute value of each row of a matrix.
"""
function maxrows(A)
	return y = vec(maximum(abs, A; dims = 2))
end
