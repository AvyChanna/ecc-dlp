from itertools import count
from multiprocessing import Process
from multiprocessing import Queue as mpQueue
from random import randrange

from .number_theory import gcd, ilog, ispower, isprime, isqrt, primegen


def pollardrho_brent(n):
	if isprime(n):
		return n
	g = n
	while g == n:
		y, c, m, g, r, q = randrange(1, n), randrange(1, n), randrange(1, n), 1, 1, 1
		while g == 1:
			x, k = y, 0
			for _ in range(r):
				y = (y**2 + c) % n
			while k < r and g == 1:
				ys = y
				for _ in range(min(m, r - k)):
					y = (y**2 + c) % n
					q = q * abs(x - y) % n
				g, k = gcd(q, n), k + m
			r *= 2
		if g == n:
			while True:
				ys = (ys**2 + c) % n
				g = gcd(x - ys, n)
				if g > 1:
					break
	return g


def pollard_pm1(n, B1=100, B2=1000):
	if isprime(n):
		return n
	m = ispower(n)
	if m:
		return m[0]
	while True:
		pg = primegen()
		q = 2
		p = next(pg)
		while p <= B1:
			q, p = pow(q, p**ilog(B1, p), n), next(pg)
		g = gcd(q - 1, n)
		if 1 < g < n:
			return g
		while p <= B2:
			q, p = pow(q, p, n), next(pg)
		g = gcd(q - 1, n)
		if 1 < g < n:
			return g
		# These bounds failed.  Increase and try again.
		B1 *= 10
		B2 *= 10


def mlucas(v, a, n):
	v1, v2 = v, (v**2 - 2) % n
	for bit in bin(a)[3:]:
		v1, v2 = ((v1**2 - 2) % n, (v1 * v2 - v) % n) if bit == "0" else ((v1 * v2 - v) % n, (v2**2 - 2) % n)
	return v1


def williams_pp1(n):
	if isprime(n):
		return n
	m = ispower(n)
	if m:
		return m[0]
	for v in count(3):
		for p in primegen():
			e = ilog(isqrt(n), p)
			if e == 0:
				break
			for _ in range(e):
				v = mlucas(v, p, n)
			g = gcd(v - 2, n)
			if 1 < g < n:
				return g
			if g == n:
				break


def multifactor(n, methods):
	def factory(method, n, output):
		output.put(method(n))

	factors = mpQueue()
	procs = [Process(target=factory, args=(m, n, factors)) for m in methods]
	for p in procs:
		p.start()
	f = factors.get()
	for p in procs:
		p.terminate()
	return f


def primefac(n, trial=1000, rho=42000, methods=(pollardrho_brent, )):
	if n < 0:
		yield -1
		n *= -1
	if n < 2:
		return
	if isprime(n):
		yield n
		return

	# Trial division
	factors, nroot = [], isqrt(n)
	for p in primegen(max(4, trial)):
		if n % p == 0:
			while n % p == 0:
				yield p
				n //= p
			nroot = isqrt(n)
		if p > nroot:
			if n != 1:
				yield n
			return
	if isprime(n):
		yield n
		return

	# Pollard's rho
	factors, difficult = [n], []
	while len(factors) != 0:
		rhocount = 0
		n = factors.pop()
		try:
			g = n
			while g == n:
				x, c, g = randrange(1, n), randrange(1, n), 1
				y = x
				while g == 1:
					if rhocount >= rho:
						raise ValueError
					rhocount += 1
					x = (x**2 + c) % n
					y = (y**2 + c) % n
					y = (y**2 + c) % n
					g = gcd(x - y, n)
			if isprime(g):
				yield g
			else:
				factors.append(g)
			n //= g
			if isprime(n):
				yield n
			else:
				factors.append(n)
		except ValueError:
			difficult.append(n)

	factors = difficult
	while len(factors) != 0:
		n = factors.pop()
		f = methods[0](n) if len(methods) == 1 else multifactor(n, methods=methods)
		if isprime(f):
			yield f
		else:
			factors.append(f)
		n //= f
		if isprime(n):
			yield n
		else:
			factors.append(n)


def factorint(n, trial=1000, rho=42000, methods=(pollardrho_brent, )):
	fac = {}
	for p in primefac(n, trial=trial, rho=rho, methods=methods):
		fac[p] = fac.get(p, 0) + 1
	return fac
