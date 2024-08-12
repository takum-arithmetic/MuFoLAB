# See LICENSE file for copyright and license details.
using SparseArrays

push!(LOAD_PATH, "src/")
import QR
import Utilities

function solve_qr(A::AbstractMatrix, b::AbstractVector)
	Q, R = QR.qr_givens(A)
	Q_full = Q * spdiagm(ones(size(A, 1)))
	z = Q_full' * b
	return R \ z
end

Utilities.run_solver_experiment(; solver = solve_qr, preconditioner = nothing)
