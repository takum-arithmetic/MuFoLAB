# See LICENSE file for copyright and license details.
module Crutches

using BFloat16s
using SoftPosit

# rounding and truncation are incomplete in SoftPosit
Base.round(p::AbstractPosit, r::RoundingMode{:ToZero})  = typeof(p)(Base.trunc(Float64(p)))
Base.round(p::AbstractPosit, r::RoundingMode{:Down})    = typeof(p)(Base.floor(Float64(p)))
Base.round(p::AbstractPosit, r::RoundingMode{:Up})      = typeof(p)(Base.ceil(Float64(p)))
Base.round(p::AbstractPosit, r::RoundingMode{:Nearest}) = typeof(p)(Base.round(Float64(p)))

Base.trunc(p::AbstractPosit) = Base.signbit(p) ? Base.ceil(p) : Base.floor(p)
Base.trunc(::Type{T}, p::AbstractPosit) where {T <: Integer} = Base.trunc(T, Float64(p))

# integer powers are broken for negative exponents, thus we hijack
# the original function and, if the exponent is negative, negate it and
# invert the result.
function Base.:(^)(x::AbstractPosit, y::Integer)
	return if sign(y) < 0
		inv(invoke(Base.:(^), Tuple{Real, Integer}, x, -y))
	else
		invoke(Base.:(^), Tuple{Real, Integer}, x, y)
	end
end

# also fix posit type promotion
Base.promote_rule(::Type{Float16}, ::Type{<:AbstractPosit}) = Float16
Base.promote_rule(::Type{Float32}, ::Type{<:AbstractPosit}) = Float32
Base.promote_rule(::Type{Float64}, ::Type{<:AbstractPosit}) = Float64

# BFloat16s defines truncation to Int, whereas it should be more generally
# defined to be Integer
Base.trunc(::Type{Integer}, f::BFloat16) = Base.trunc(Int, f)

end
