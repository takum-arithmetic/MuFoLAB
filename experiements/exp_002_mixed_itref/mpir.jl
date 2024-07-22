"""
    Experiement 2: Mixed Iterative Refinement (no scaling/preconditioning)
    Solve Ax = b using mixed-precision iterative refinement (MPIR) with
    different types.  See Experiment 1 for comparison of solutions.  

    Overview: perform computationally expensive operations in low
    precision and high precision to compute the residual to obtain
    a more accurate solution [1].

    MPIR: Three precisions are used: high, medium, and low.  The high
    is used to compute the residual, the medium is working precision.  
    Low precision is used for the O(n^3) operations (i.e., LU 
    factorization).

    Reference:
    [1] Forsythe, G. E., & Moler, C. B. (1967). Computer solution of linear 
    algebraic systems.
"""

include("configs.jl")

# Run the experiment
# 1. Get matrix
# 2. Solve Ax = b using MPIR
# 3. Display results
