function sparsecast(t::Type{T}, M::AbstractMatrix) where {T <: AbstractFloat}
	if typeof(M) <: SparseMatrixCSC
		return SparseMatrixCSC{t, Int64}(M)
	elseif typeof(M) <: Symmetric
		return Symmetric{t, SparseMatrixCSC{t, Int64}}(M)
	else
		throw("Unhandled sparse matrix type")
	end
end

function xgesol(i::Int64, t::Type{T}) where {T <: AbstractFloat}
	A = sparsecast(t, matrixdepot(sp(i)))
	y = ones(t, size(A, 1))
	b = A * y
	try
		x = A \ b
		return norm(Float64.(y) - Float64.(x), Inf)
	catch
		return NaN
	end
end
