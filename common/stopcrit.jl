"""
    STOPCRIT(A, b, x, tol) is a criterion to stop the iterative method.
    The criterion is satisfied if [1, p. 54]:

    ||r|| ≤ tol ⋅ ( ||b|| + ||A|| ⋅ ||x|| )

    This is equivalent to backward errors satistfying:
     ||ΔA|| ≤ tol * ||A||  and ||Δb|| ≤ tol * ||b||. 

    The criterion yields foward error bound:
    ||e|| ≤ ||A^{-1}|| ⋅ ||r|| ≤ tol ⋅ ||A^{-1}|| ⋅ (||A||⋅||x|| + ||b||)

    Reference:
    [1] Barrett, Richard, et al. Templates for the solution of linear 
    systems: building blocks for iterative methods. Society for Industrial 
    and Applied Mathematics, 1994.

    @param A: Matrix
    @param b: Vector
    @param x: Vector
    @param tol: tolerance (default: 1e-06)
    @return true if the criterion is satisfied, false otherwise

    Example:
    ```julia-repl
    julia> A = [1.0 2.0; 3.0 4.0]
    julia> b = [1.0, 2.0]
    julia> x = [1.0, 1.0]
    julia> stopcrit(A, b, x, 1e-06)
    false
    ```
"""
function stopcrit(A, b, x, tol = 1e-06)
	T = typeof(A[1, 1])
	r = b - A * x
	N = infnorm
	return N(r) <= T.(tol) * (N(b) + N(A) * N(x)) ? true : false
end
