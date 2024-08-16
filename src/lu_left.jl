# See LICENSE file for copyright and license details.
using LinearAlgebra
using SparseArrays

function lu_left(A::SparseMatrixCSC{T, Int64}) where {T <: AbstractFloat}
    n = size(A, 1)
    P = sparse(I,n,n)
    L = spzeros(T,n,n)
    U = spzeros(T,n,n)
    for k = 1:n
        x = [ L[:,1:k-1] [ zeros(k-1,n-k+1) ; eye(n-k+1) ]] \ (P * A[:,k])
        U[1:k-1,k] = x[1:k-1]
        a = maximum(abs.(x[k:n]))
        i = argmax(abs.(x[k:n]))
        i = i + k - 1
        L[[i k],:] = L[[k i], :]
        P[[i k],:] = P[[k i], :]
        x[[i k]] = x[[k i]]
        U[k,k] = x[k]
        L[k,k] = 1
        L[k+1:n,k] = x[k+1:n] / x[k]
    end
    return L, U, P
end

function eye(n)
    return sparse(I,n,n)
end