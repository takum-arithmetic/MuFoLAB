# See LICENSE file for copyright and license details.
module Crutches

using BFloat16s
using LinearAlgebra
using MicroFloatingPoints

# BFloat16s defines truncation to Int, whereas it should be more generally
# defined to be Integer
Base.trunc(::Type{Integer}, f::BFloat16) = Base.trunc(Int, f)

# define integer conversion for microfloats
Integer(x::Floatmu{szE, szf}) where {szE, szf} = Integer(Float32(x))

# overwrite nameof for microfloats
function Base.nameof(::Type{Floatmu{szE, szf}}) where {szE, szf}
	return "Float" * string(1 + szE + szf) * "_" * string(szE)
end

end
