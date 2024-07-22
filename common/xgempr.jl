function xgempr(
        A::AbstractMatrix{T}, 
        b::Vector{T}, 
        t::Vector{DataType}, 
        maxits::Int = 20
    ) where T <: AbstractFloat
    
    niters = 1
    x = zeros(T, lastindex(b))
    
    F = lu(t[1].(A), check = true)
    L = convert(Matrix{T}, F.L)
    U = convert(Matrix{T}, F.U)
    P = convert(Matrix{T}, F.P)
  
    while stopcrit(A, b, x) == false && niters < maxits
        niters += 1
        r = t[3].(b) - t[3].(A) * t[3].(x)
        x = x + U \ (L \ (P * T.(r))) 
    end
    return x, niters
end