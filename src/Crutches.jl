# See LICENSE file for copyright and license details.
module Crutches

using BFloat16s
using Float8s
using LinearAlgebra

# BFloat16s defines truncation to Int, whereas it should be more generally
# defined to be Integer
Base.trunc(::Type{Integer}, f::BFloat16) = Base.trunc(Int, f)

# Float8s does not define eps()
Base.eps(::Type{Float8}) = Float8(2^-4)

# Integer truncation is incomplete in Float8s
Base.trunc(::Type{T}, f::Float8) where {T <: Integer} = Base.trunc(T, Float64(f))

# givens_algorithm is broken in Float8s given its small dynamic range,
# we add the second count >= 20 check while pursuing an upstream fix.
#
# The function is licensed as follows:
#
# Copyright (c) 2009-2024: Jeff Bezanson, Stefan Karpinski, Viral B. Shah,
# and other contributors: https://github.com/JuliaLang/julia/contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
function LinearAlgebra.givensAlgorithm(f::Float8, g::Float8)
	T = Float8
	onepar = one(T)
	T0 = typeof(onepar) # dimensionless
	zeropar = T0(zero(T)) # must be dimensionless

	# need both dimensionful and dimensionless versions of these:
	safmn2 = LinearAlgebra.floatmin2(T0)
	safmn2u = LinearAlgebra.floatmin2(T)
	safmx2 = one(T) / safmn2
	safmx2u = oneunit(T) / safmn2

	if g == 0
		cs = onepar
		sn = zeropar
		r = f
	elseif f == 0
		cs = zeropar
		sn = onepar
		r = g
	else
		f1 = f
		g1 = g
		scalepar = max(abs(f1), abs(g1))
		if scalepar >= safmx2u
			count = 0
			while true
				count += 1
				f1 *= safmn2
				g1 *= safmn2
				scalepar = max(abs(f1), abs(g1))
				if scalepar < safmx2u || count >= 20
					break
				end
			end
			r = sqrt(f1 * f1 + g1 * g1)
			cs = f1 / r
			sn = g1 / r
			for i in 1:count
				r *= safmx2
			end
		elseif scalepar <= safmn2u
			count = 0
			while true
				count += 1
				f1 *= safmx2
				g1 *= safmx2
				scalepar = max(abs(f1), abs(g1))
				if scalepar > safmn2u || count >= 20
					break
				end
			end
			r = sqrt(f1 * f1 + g1 * g1)
			cs = f1 / r
			sn = g1 / r
			for i in 1:count
				r *= safmn2
			end
		else
			r = sqrt(f1 * f1 + g1 * g1)
			cs = f1 / r
			sn = g1 / r
		end
		if abs(f) > abs(g) && cs < 0
			cs = -cs
			sn = -sn
			r = -r
		end
	end
	return cs, sn, r
end

end
