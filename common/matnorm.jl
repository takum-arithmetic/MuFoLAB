"""
    matnorm(A,p=Inf) = ||A||_p

    @param A: Matrix
    @param p: norm type, p = 1 or p = Inf (default: Inf) 
    @return max_{x!=0} ||Ax||_p/||x||_p
"""
function matnorm(A::AbstractMatrix{T}, p = Inf) where {T<:AbstractFloat}
	v = p == 1 ? maximum(sum(abs.(A), dims = 1)) : maximum(sum(abs.(A), dims = 2))
	return v
end
