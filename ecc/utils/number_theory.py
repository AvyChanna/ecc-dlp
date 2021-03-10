from itertools import compress, takewhile
from math import log


def gcd(a, *r):
	if not isinstance(a, int):
		r = iter(a)
		a = next(r)
	for b in r:
		while b:
			a, b = b, a % b
	return abs(a)


def xgcd(a, b):
	if a == 0 and b == 0: return (0, 0, 1)
	if a == 0: return (abs(b), 0, b // abs(b))
	if b == 0: return (abs(a), a // abs(a), 0)
	x_sign = 1
	y_sign = 1
	if a < 0:
		a = -a
		x_sign = -1
	if b < 0:
		b = -b
		y_sign = -1
	x, y, r, s = 1, 0, 0, 1
	while b != 0:
		q = a // b
		a, b, r, s, x, y = b, a % b, x - q * r, y - q * s, r, s
	return (a, x * x_sign, y * y_sign)


def modinv(a, m):
	g, x, _ = xgcd(a, m)
	if g != 1: return None
	return x % m


def iterprod(l):
	z = 1
	for x in l:
		z *= x
	return z


def crt(rems, mods):
	if len(mods) == 1: return rems[0]
	N = iterprod(mods)
	return sum(r * (N // m) * modinv(N // m, m) for (r, m) in zip(rems, mods) if m != 1) % N


def legendre(a, p):
	return ((pow(a, (p - 1) >> 1, p) + 1) % p) - 1


def jacobi(a, n):
	if (n % 2 == 0) or (n < 0): return None  # n must be a positive odd number     TODO delete this check?
	if (a == 0) or (a == 1): return a
	a, t = a % n, 1
	while a != 0:
		while not a & 1:
			a //= 2
			if n & 7 in (3, 5): t *= -1
		a, n = n, a
		if (a & 3 == 3) and (n & 3) == 3: t *= -1
		a %= n
	return t if n == 1 else 0


def kronecker(a, n):
	if n == -1: return -1 if a < 0 else 1
	if n == 0: return 1 if abs(a) == 1 else 0
	if n == 1: return 1
	if n == 2: return 0 if a % 2 == 0 else (1 if a % 8 in [1, 7] else -1)
	if n < 0: return kronecker(a, -1) * kronecker(a, -n)
	f = 0
	while n % 2 == 0:
		n //= 2
		f += 1
	return kronecker(a, 2)**f * jacobi(a, n)


def sqrtmod_prime(a, p):
	a %= p
	if p % 4 == 3: return pow(a, (p + 1) >> 2, p)
	elif p % 8 == 5:
		v = pow(a << 1, (p - 5) >> 3, p)
		return (a * v * (((a * v * v << 1) % p) - 1)) % p
	elif p % 8 == 1:
		if a == 0:
			return 0
		h = 2
		while legendre(h * h - 4 * a, p) != -1:
			h += 1
		k, v, w, q, Q = (p + 1) // 2, 2, h % p, 1, 1
		for kj in bin(k)[2:]:
			q = (q * Q) % p
			if kj == '1':
				Q, v, w = (q * a) % p, (w * v - h * q) % p, (w * w - 2 * q * a) % p
			else:
				Q, w, v = q, (w * v - h * q) % p, (v * v - 2 * q) % p
		return (v * k) % p
	else:
		return a  # p = 2


def sprp(n, b):
	t, s = (n - 1) // 2, 1
	while t % 2 == 0:
		t //= 2
		s += 1
	#assert 1 + 2**s * t == n
	x = pow(b, t, n)
	if x == 1 or x == n - 1: return True
	for j in range(1, s):
		x = pow(x, 2, n)
		if x == 1: return False
		elif x == n - 1: return True
	return False


def mrab(n, basis):
	return all(sprp(n, b) for b in basis)


def miller(n):
	return n == 2 if n % 2 == 0 or n < 3 else mrab(n, range(2, min(n - 1, int(2 * log(n)**2))))


def isprime(n):
	return miller(n)


def isqrt(n):
	if n < 0: return int(n)
	c = n * 4 // 3
	d = c.bit_length()
	a = d >> 1
	if d & 1:
		x = 1 << a
		y = (x + (n >> a)) >> 1
	else:
		x = (3 << a) >> 2
		y = (x + (c >> a)) >> 1
	if x != y:
		x, y = y, (y + n // y) >> 1
		while y < x:
			x, y = y, (y + n // y) >> 1
	return x


def introot(n, r=2):  # TODO Newton iteration?
	if n < 0: return None if r % 2 == 0 else -introot(-n, r)
	if n < 2: return n
	if r == 1: return n
	if r == 2: return isqrt(n)
	#if r % 2 == 0: return introot(isqrt(n), r//2)      # TODO Check validity of this line.
	lower = upper = 1 << (n.bit_length() // r)
	while lower**r > n:
		lower >>= 2
	while upper**r <= n:
		upper <<= 2
	while lower != upper - 1:
		mid = (lower + upper) // 2
		m = mid**r
		if m == n: return mid
		elif m < n: lower = mid
		elif m > n: upper = mid
	return lower


def primegen(limit=None):
	if limit is None:
		limit = float('inf')
	yield from takewhile(lambda x: x < limit, (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47))
	pl, pg = [3, 5, 7], primegen()
	for p in pl:
		next(pg)
	n = next(pg)
	nn = n * n
	while True:
		n = next(pg)
		ll, nn = nn, n * n
		sl = (nn - ll)
		sieve = bytearray([True]) * sl
		for p in pl:
			k = (-ll) % p
			sieve[k::p] = bytearray([False]) * ((sl - k) // p + 1)
		if nn > limit: break  # TODO bring this condition up to the while statement
		yield from (ll + k for k in compress(range(0, sl, 2), sieve[::2]))
		pl.append(n)
	yield from takewhile(lambda x: x < limit, (ll + k for k in compress(range(0, sl, 2), sieve[::2])))


def ispower(n, r=0):
	if r == 0:
		if n in (0, 1, -1): return (n, 1)
		for r in primegen(n.bit_length() + 1):
			x = ispower(n, r)
			if x is not None: return (x, r)
		return None
	# TODO tricks for special cases
	if (r == 2) and (n & 2): return None
	if (r == 3) and (n & 7) in (2, 4, 6): return None
	x = introot(n, r)
	return None if x is None else (x if x**r == n else None)


def ilog(x, b):
	l = 0
	while x >= b:
		x //= b
		l += 1
	return l
