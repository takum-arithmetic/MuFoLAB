# See LICENSE file for copyright and license details.
module Float128Conversions

using BFloat16s
using Float8s
using Quadmath
using SoftPosit
using Takums

# this module facilitates conversions between the numerical types
# under test and quadruple precision floating-point numbers

# For takum it is more difficult as we cannot fall back to the
# floating-point functions, given takum64 is more precise than
# float64. Use takum64 as a common ground and specify the
# conversions from and to takum64 in a special function.
function Takums.Takum64(x::Float128)
	# catch special cases early on
	if !isfinite(x)
		return NaR64
	elseif iszero(x)
		return zero(Takum64)
	end

	# We can now assume that x has the regular form
	#
	#	x = (-1)^S * sqrt(e)^l
	#
	# and the first step is to determine s and l from x.
	S = Integer(x < 0)
	l = 2 * log(abs(x))

	# Clamp l to representable exponents
	bound = Float128(BigFloat("254.999999999999999777955395074968691915273666381835938"))
	l = (l < -bound) ? -bound : (l > bound) ? bound : l

	# It holds l = (-1)^s (c + m), where c is the characteristic
	# and m the mantissa, both quantities directly encoded in the
	# takum modulo a possible final negation.
	cpm = (S == 0) ? l : -l
	c = Integer(floor(cpm))
	m = cpm - c

	# determine D
	D = Integer(c >= 0)

	# determine r and R
	r = (D == 0) ? Unsigned(floor(log2(-c))) : Unsigned(floor(log2(c + 1)))
	R = (D == 0) ? 7 - r : r

	# determine the characteristic bits C
	C = (D == 0) ? (c + 2^(r + 1) - 1) : (c - 2^r + 1)

	# determine precision p, i.e. the number of mantissa bits
	p = 64 - 5 - r

	# We extract the lower 112 bits of m, the significand bits.
	#
	# If the exponent bits are not all zero (indicating subnormal or
	# zero), we apply the implicit 1 bit and then apply the
	# m-exponent by shifting -exponent(m) to the right
	# (exponent(m) is always negative as m in [0,1).

	# first just get the float128 bits
	M = reinterpret(UInt128, m)

	# extract the coded exponent
	coded_exponent = (M & UInt128(0x7fff_0000_0000_0000_0000_0000_0000_0000)) >> 112

	# now store only the significand bits
	M &= UInt128(0x0000_ffff_ffff_ffff_ffff_ffff_ffff_ffff)

	if (coded_exponent != 0)
		# apply the implicit 1 bit
		M |= UInt128(1) << 112

		# shift M to the left such that we have no gaps
		M <<= 15

		# shift M to the right by -exponent(m) + 1 (-1, because
		# we also shift out the implicit 1 bit so it's truly
		# a 0. representation
		M >>= -exponent(m) - 1
	else
		# no implicit 1 bit, we just directly move to the left
		M <<= 16
	end

	# shift M to the right by 2 + 3 + r to accomodate for
	# the takum sign, direction bit and exponent
	M >>= 2 + 3 + r

	# assemble 128-bit takum
	t128 =
		(UInt128(S) << 127) |
		(UInt128(D) << 126) |
		(UInt128(R) << 123) |
		(UInt128(C) << (123 - r)) |
		M

	# round to 64-bit, adhering to proper saturation
	t64 = reinterpret(
		Takum64,
		UInt64(t128 >> 64) + UInt64((t128 & (UInt128(1) << 63)) >> 63),
	)

	if (iszero(t64) && !iszero(x))
		if x < 0
			# overflow to 0
			t64 = reinterpret(Takum64, Int64(-1))
		else
			# underflow to 0
			t64 = reinterpret(Takum64, Int64(1))
		end
	elseif (isnan(t64))
		if x < 0
			# underflow to NaR
			t64 = reinterpret(Takum64, typemin(Int64) + 1)
		else
			# overflow to NaR
			t64 = reinterpret(Takum64, typemax(Int64))
		end
	end

	return t64
end

Takums.Takum8(x::Float128) = Takum8(Takum64(x))
Takums.Takum16(x::Float128) = Takum16(Takum64(x))
Takums.Takum32(x::Float128) = Takum32(Takum64(x))

function Quadmath.Float128(t::Takum64)
	# reinterpret the takum as an unsigned 64-bit integer
	T = reinterpret(UInt64, t)

	# get the obvious bits
	S = (T & (UInt64(1) << 63)) != 0
	D = (T & (UInt64(1) << 62)) != 0
	R = (T & (UInt64(7) << 59)) >> 59
	r = (D == 0) ? (7 - R) : R

	# shift to the left, shift to the right, obtain C and M without
	# bitmasks
	C = (T << 5) >> (64 - r)
	M = T << (5 + r)

	# obtain c from C
	c = if (D == 0)
		(-Int16(2)^(r + 1) + Int16(1) + Int16(C))
	else
		(Int16(2)^r - Int16(1) + Int16(C))
	end

	# build a fixed point representation of (-1)^S * l = c + m
	cM = Int128(M) | (Int128(c) << 64)

	# this representation has at most 64+9 = 73 significant digits,
	# way below the limits of what float128 can represent.
	# Cast to float128 and divide by 2^64 to obtain c + m
	cpm = Float128(cM) / Float128(2.0)^64

	# l follows directly
	l = (S == 0) ? cpm : -cpm

	# determine sqrt(e)^l
	lraised = exp(0.5 * l)

	return (S == 0) ? lraised : -lraised
end

Quadmath.Float128(t::Takum8) = Float128(Takum64(t))
Quadmath.Float128(t::Takum16) = Float128(Takum64(t))
Quadmath.Float128(t::Takum32) = Float128(Takum64(t))

# For posit use Float64 as a common ground, as we cannot go
# higher than posit32 anyway
SoftPosit.Posit8(x::Float128) = Posit8(Float64(x))
SoftPosit.Posit16(x::Float128) = Posit16(Float64(x))
SoftPosit.Posit32(x::Float128) = Posit32(Float64(x))
Quadmath.Float128(x::Posit8) = Float128(Float64(x))
Quadmath.Float128(x::Posit16) = Float128(Float64(x))
Quadmath.Float128(x::Posit32) = Float128(Float64(x))

# For bfloat16 use Float32 as a common ground
BFloat16s.BFloat16(x::Float128) = BFloat16(Float32(x))
Quadmath.Float128(x::BFloat16) = Float128(Float32(x))

# For float8 use Float32 as a common ground, as well
Float8s.Float8(x::Float128) = Float8(Float32(x))
Quadmath.Float128(x::Float8) = Float128(Float32(x))

end
