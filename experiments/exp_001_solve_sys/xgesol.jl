function xgesol(i::Int64, t::Type{T}) where {T <: AbstractFloat}
	A = Array(t.(matrixdepot(sp(i))))
	y = t.(ones(size(A, 1)))
	b = A * y
	try
		x = A \ b
		return round(Float64.(infnorm(y - x)); digits = 6)
	catch
		return NaN
	end
end
