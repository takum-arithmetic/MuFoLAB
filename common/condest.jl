""" 
    condest(A::AbstractMatrix{T}) where {T <: AbstractFloat} -> Float64 

    Estimates the condition number of a matrix `A` using the p-norm.
    condest(A) = ||A||_p * ||A^{-1}||_p.  
    Method based on [1, p.372].  Note: [1] uses the 1-norm.

Reference 
----------
[1] Cline, A. K., Moler, C. B., Stewart, G. W., & Wilkinson, J. H. (1979). 
An estimate for the condition number of a matrix. SIAM Journal on 
Numerical Analysis, 16(2), 368-375.
"""
function condest(A::AbstractMatrix{T}; p = Inf) where {T <: AbstractFloat}
        n = size(A, 1)
        b = T.(rand([1, -1], n))
        x = A' \ b
        y = A \ x
        return norm(A, p) * (norm(y, p) / norm(x, p))
end
