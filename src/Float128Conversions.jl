# See LICENSE file for copyright and license details.
module Float128Conversions

using BFloat16s
using Quadmath
using MicroFloatingPoints
using Posits
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
		return NaRTakum64
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
	# catch special cases
	if isnan(t)
		return Float128(NaN)
	elseif t == zero(Takum64)
		return zero(Float128)
	end

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

# do similarly for linear takums
function Takums.LinearTakum64(x::Float128)
	# catch special cases early on
	if !isfinite(x)
		return NaRLinearTakum64
	elseif iszero(x)
		return zero(LinearTakum64)
	end

	# We can now assume that x has the regular form
	#
	#	x = (-1)^S * (1+g) * 2^h
	#	x = [(1-3S) + f] * 2^e, e = (-1)^S (c + S)

	# get g and h, but normalise as frexp() returns g in [0.5,1),
	# not [1,2), and we actually want the form (1+g) where g
	# in [0,1)
	g, h = frexp(abs(x))
	g *= 2
	h -= 1
	g -= 1

	# and the first step is to determine c and m from x, so we just
	# treat it as if it was a normal logarithmic takum and carry
	# on from c and m
	S = Integer(x < 0)

	if S == 0
		c = h
		m = g
	else
		if g == 0
			c = -h
			m = Float128(0)
		else
			c = -h - 1
			m = 1 - g
		end
	end

	# clamp c
	c = (c < -255) ? -255 : (c > 254) ? 254 : c

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

	# assemble 128-bit linear takum
	t128 =
		(UInt128(S) << 127) |
		(UInt128(D) << 126) |
		(UInt128(R) << 123) |
		(UInt128(C) << (123 - r)) |
		M

	# round to 64-bit, adhering to proper saturation
	t64 = reinterpret(
		LinearTakum64,
		UInt64(t128 >> 64) + UInt64((t128 & (UInt128(1) << 63)) >> 63),
	)

	if (iszero(t64) && !iszero(x))
		if x < 0
			# overflow to 0
			t64 = reinterpret(LinearTakum64, Int64(-1))
		else
			# underflow to 0
			t64 = reinterpret(LinearTakum64, Int64(1))
		end
	elseif (isnan(t64))
		if x < 0
			# underflow to NaR
			t64 = reinterpret(LinearTakum64, typemin(Int64) + 1)
		else
			# overflow to NaR
			t64 = reinterpret(LinearTakum64, typemax(Int64))
		end
	end

	return t64
end

Takums.LinearTakum8(x::Float128) = LinearTakum8(LinearTakum64(x))
Takums.LinearTakum16(x::Float128) = LinearTakum16(LinearTakum64(x))
Takums.LinearTakum32(x::Float128) = LinearTakum32(LinearTakum64(x))

function Quadmath.Float128(t::LinearTakum64)
	# catch special cases
	if isnan(t)
		return Float128(NaN)
	elseif t == zero(LinearTakum64)
		return zero(Float128)
	end

	# reinterpret the linear takum as an unsigned 64-bit integer
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

	# obtain f (which is equivalent to m) from M
	f = Float128(M) / Float128(2.0)^64

	# compute e
	if S == 0
		e = c
	else
		e = -(c + 1)
	end

	# compute the linear takum value and return
	return (Float128(1 - 3 * S) + f) * Float128(2.0)^e
end

Quadmath.Float128(t::LinearTakum8) = Float128(LinearTakum64(t))
Quadmath.Float128(t::LinearTakum16) = Float128(LinearTakum64(t))
Quadmath.Float128(t::LinearTakum32) = Float128(LinearTakum64(t))

# posit
function Posits.Posit64(x::Float128)
	# catch special cases early on
	if !isfinite(x)
		return NaRPosit64
	elseif iszero(x)
		return zero(Posit64)
	end

	# We can now assume that x has the regular form
	#
	#	x = (-1)^S * (1+g) * 2^h
	#	x = [(1-3S) + f] * 2^e, e = (-1)^S (c + S)

	# get g and h, but normalise as frexp() returns g in [0.5,1),
	# not [1,2), and we actually want the form (1+g) where g
	# in [0,1)
	g, h = frexp(abs(x))
	g *= 2
	h -= 1
	g -= 1

	# and the first step is to determine c and m from x
	S = Integer(x < 0)

	if S == 0
		c = h
		m = g
	else
		if g == 0
			c = -h
			m = Float128(0)
		else
			c = -h - 1
			m = 1 - g
		end
	end

	# clamp c (arbitrary, but I don't care. It fits posit64, that's enough)
	c = (c < -255) ? -255 : (c > 254) ? 254 : c

	# c = 4r+e, r = (R0 = 0) ? -k : k-1
	R0 = (c >= 0)
	R0b = !R0
	k = 0

	if R0 == 0
		# c is negative, c = -4k + e <-> -c = 4k - e -> -c+3 = 4k + (3-e)
		tmp = -c + 3
		e = 3 - (tmp % 4)
		k = Int(floor(tmp / 4))
	else
		# c is positive, c = 4(k-1) + e -> c = 4k + e - 4 -> c + 4 = 4k + e
		tmp = c + 4
		e = tmp % 4
		k = Int(floor(tmp / 4))
	end

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

	# get precision
	p = 128 - (1 + k + 1 + 2)

	# shift M to the right by (n - p) to accomodate for
	# the posit sign, regime and exponent bits
	M >>= (128 - p)

	if R0 == 1
		R = UInt128(2.0)^k - 1
	else
		R = 0
	end

	# assemble 128-bit posit
	t128 =
		(UInt128(S) << 127) |
		(R << (p + 3)) |
		(UInt128(R0b) << (p + 2)) |
		(UInt128(e) << p) |
		M

	# round to 64-bit, adhering to proper saturation
	t64 = reinterpret(
		Posit64,
		UInt64(t128 >> 64) + UInt64((t128 & (UInt128(1) << 63)) >> 63),
	)

	if (iszero(t64) && !iszero(x))
		if x < 0
			# overflow to 0
			t64 = reinterpret(Posit64, Int64(-1))
		else
			# underflow to 0
			t64 = reinterpret(Posit64, Int64(1))
		end
	elseif (isnan(t64))
		if x < 0
			# underflow to NaR
			t64 = reinterpret(Posit64, typemin(Int64) + 1)
		else
			# overflow to NaR
			t64 = reinterpret(Posit64, typemax(Int64))
		end
	end

	return t64
end

Posits.Posit8(x::Float128) = Posit8(Float64(x))
Posits.Posit16(x::Float128) = Posit16(Float64(x))
Posits.Posit32(x::Float128) = Posit32(Float64(x))

function Quadmath.Float128(t::Posit64)
	# catch special cases
	if isnan(t)
		return Float128(NaN)
	elseif t == zero(Posit64)
		return zero(Float128)
	end

	# reinterpret the posit as an unsigned 64-bit integer
	T = reinterpret(UInt64, t)

	# get the sign
	S = (T & (UInt64(1) << 63)) != 0

	# Posits have the form
	#
	#      SR...RNEEM...M
	#
	# where N = NOT(R) and there are k R's, with k >= 1. If we shift the
	# posit one to the left we obtain
	#
	#      R...RNEEM...M0
	#
	# and it holds, when we XOR this with it once more shifted to the left
	#
	#     R...RNEEM...M0
	#     R..RNEEM...M00 XOR
	#     ------------------
	#     00001xxxxxxxx0
	#
	# If we count the leading zeros of this result we directly obtain k-1.
	#
	# By handling the special cases beforehand we know that the result of
	# (p << 1) ^ (p << 2) is never zero.
	k = leading_zeros((T << 1) ⊻ (T << 2)) + 1

	# Determine c = 4r+e with r = (R0 = 0) ? -k : k-1
	R0 = ((T & UInt64(0x4000000000000000)) != 0) ? 1 : 0

	# Determine exponent value
	if k <= 64 - 4
		# both exponent bits are explicitly given
		e = (T >> (64 - 4 - k)) & UInt64(0x3)
	elseif k == 64 - 4 + 1
		# only the higher exponent bit is given, right at the end
		e = 2 * (T & UInt64(0x1))
	else
		# no exponent bit is given, ghost bits imply zero
		e = 0
	end

	c = ((R0 == 0) ? (-4 * Int(k)) : (4 * (Int(k) - 1))) + Int(e)

	# the mantissa bits follow by shifting p by k+4 to the left
	M = (k < 64 - 4) ? (T << (k + 4)) : 0

	# obtain f (which is equivalent to m) from M
	f = Float128(M) / Float128(2.0)^64

	# compute exponent
	if S == 0
		exponent = c
	else
		exponent = -(c + 1)
	end

	# compute the linear takum value and return
	return (Float128(1 - 3 * S) + f) * Float128(2.0)^exponent
end

Quadmath.Float128(x::Posit8) = Float128(Float64(x))
Quadmath.Float128(x::Posit16) = Float128(Float64(x))
Quadmath.Float128(x::Posit32) = Float128(Float64(x))

# For bfloat16 use Float32 as a common ground
BFloat16s.BFloat16(x::Float128) = BFloat16(Float32(x))
Quadmath.Float128(x::BFloat16) = Float128(Float32(x))

# For microfloats use Float32 as a common ground, as well
Quadmath.Float128(x::Floatmu{szE, szf}) where {szE, szf} = Float128(Float32(x))
MicroFloatingPoints.Floatmu{szE, szf}(x::Float128) where {szE, szf} = Floatmu{szE, szf}(Float32(x))

end
