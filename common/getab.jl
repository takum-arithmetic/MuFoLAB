function getab(i::Int16, t::Type{T}) where {T<:AbstractFloat}
    A = Matrix(t.(matrixdepot(sp(i))))
    b = A * t.(ones(size(A,1)))
    return A, b
end