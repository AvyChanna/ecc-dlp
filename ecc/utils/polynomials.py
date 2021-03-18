from itertools import zip_longest
from random import randrange

from .number_theory import isprime, legendre, modinv, sqrtmod_prime


def polyaddmodp(a, b, p):
	if a is None or b is None:
		return None
	c = [(x + y) % p for (x, y) in zip_longest(a, b, fillvalue=0)]
	while len(c) > 0 and c[-1] == 0:
		del c[-1]
	return c


def polysubmodp(a, b, p):
	if a is None or b is None:
		return None
	c = [(x - y) % p for (x, y) in zip_longest(a, b, fillvalue=0)]
	while len(c) > 0 and c[-1] == 0:
		del c[-1]
	return c


def polymulmodp(a, b, p):
	if a is None or b is None:
		return None
	c = [0] * (len(a) + len(b) - 1)
	for (k, x) in enumerate(a):
		for (l, y) in enumerate(b):
			c[k + l] += x * y
	for (k, x) in enumerate(c):
		c[k] = x % p
	while len(c) > 0 and c[-1] == 0:
		del c[-1]
	return c


def polydivmodmodp(a, b, p):
	if a is None or b is None or b == []:
		return (None, None)
	ndeg, ddeg = len(a) - 1, len(b) - 1
	if ndeg < ddeg:
		return ([], a[:])
	num, den = a[:], b[:]
	c = modinv(b[-1], p)
	if c is None:
		return (None, None)
	term = [0] * (ndeg - ddeg) + [(num[-1] * c) % p]
	newnum = polysubmodp(num, polymulmodp(den, term, p), p)
	quo, rem = polydivmodmodp(newnum, den, p)
	quo = polyaddmodp(term, quo, p)
	while len(quo) > 0 and quo[-1] == 0:
		del quo[-1]
	while len(rem) > 0 and rem[-1] == 0:
		del rem[-1]
	return (quo, rem)


def polypowmodpmodpoly(a, e, p, f):  # a**e mod p mod f
	if a is None or f is None or f == []:
		return None
	ans, a = [1], polydivmodmodp(a, f, p)[1]
	for k in bin(e)[2:]:
		ans = polydivmodmodp(polymulmodp(ans, ans, p), f, p)[1]
		if k == '1':
			ans = polydivmodmodp(polymulmodp(ans, a, p), f, p)[1]
	return ans


def gcmd(f, g, p):
	if (f is None) or (g is None) or f == [] == g:
		return None
	df, dg = len(f) - 1, len(g) - 1
	if dg > df:
		u, v = g[:], f[:]
	else:
		u, v = f[:], g[:]
	while v != []:
		r = polydivmodmodp(u, v, p)[1]
		if r is None:
			return None
		u, v = v[:], r
	c = modinv(u[-1], p)
	if c is None:
		return None
	return [(x * c) % p for x in u]


def polyval(f, x, m=None):
	out = 0
	if m is None:
		for a in reversed(f):
			out = a + out * x
	else:
		for a in reversed(f):
			out = (a + out * x) % m
	return out


def polyroots_prime(f, p, sqfr=False):
	assert isprime(p)
	g = f[:]
	while len(g) > 0 and g[-1] % p == 0:
		del g[-1]
	if g == []:
		yield from range(p)
		return
	if not sqfr:
		g = gcmd(g, polysubmodp(polypowmodpmodpoly([0, 1], p, p, g), [0, 1], p), p)
	else:
		g = [x % p for x in g]
		while len(g) > 0 and g[-1] == 0:
			del g[-1]
	assert not (g is None), (g, p)
	if g == []:
		yield from range(p)
		return
	if g[0] == 0:
		yield 0
		while g[0] == 0:
			del g[0]
		if g == []:
			assert False
	if len(g) == 1:
		return
	if len(g) == 2:
		yield (-g[0] * modinv(g[1], p)) % p
		return  # Linear.
	if len(g) == 3:  # Quadratic.
		if p == 2:
			yield from (x for x in (0, 1) if polyval(g, x, 2) == 0)
			return
		c, b, a = g  # ax^2  +  bx  +  c
		ai = modinv(a, p)
		c, b = (c * ai) % p, (b * ai) % p
		d = b * b - 4 * c
		l = legendre(d, p)
		if l == -1:
			return
		sq = sqrtmod_prime(d, p)
		inv2 = (p + 1) // 2  # inv2 == modinv(2, p)
		yield from {((-b + sq) * inv2) % p, ((-b - sq) * inv2) % p}
		return
	h = [1]
	while len(h) in (1, len(g)):
		h = gcmd(polysubmodp(polypowmodpmodpoly([randrange(p), 1], (p - 1) // 2, p, g), [1], p), g, p)
	q, s = polydivmodmodp(g, h, p)
	assert s == []
	yield from polyroots_prime(h, p, True)
	yield from polyroots_prime(q, p, True)
	return
