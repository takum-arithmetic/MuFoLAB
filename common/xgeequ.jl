"""
    Equilibrate a matrix `A` by scaling rows and columns to have unit norms.

    @param A: the matrix to equilibrate
    @return: the equilibrated matrix

    Note: A not mutable since custom types defined without the mutable 
    keyword are not mutable by default.
"""
function xgeequ!(A::AbstractMatrix{T}) where {T}
	R = 1 ./ maxrows(A)
	A = R .* A
	S = 1 ./ maxcols(A)'
	A = A .* S
	return A, R
end
