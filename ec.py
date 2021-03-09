from collections import namedtuple
from textwrap import dedent

from utils import isprime, modinv, polyroots_prime, sqrtmod_prime

ec_point = namedtuple("ec_point", "x y z", defaults=(1, ))
ec_point.is_origin = (lambda x: x.z == 0)  # type: ignore
ec_point.__repr__ = (lambda x: f"Pt({x.x}, {x.y})" if x.z == 1 else "Pt(Origin)")  # type: ignore


class ec_curve:
	# Elliptic curve over finite prime field Z/nZ
	def __init__(self, a: int, b: int, p: int):
		assert isprime(p)
		self.p = p
		self.a = a % p
		self.b = b % p
		a4 = 4 * (self.a**3)
		b27 = 27 * (self.b**2)
		a4b27 = a4 + b27
		self.discriminant = -16 * a4b27
		assert self.discriminant != 0
		self.j_invariant = 1728 + (a4 * modinv(a4b27, self.p))
		self.card = None

	def from_x(self, x):
		# y2 = [x^3 + ax +b] mod p
		y2 = (x * x * x + self.a * x + self.b) % self.p
		# screw tonelli shanks, hensel lift rocks
		candidates = [sqrtmod_prime(y2, self.p)]
		assert len(candidates) != 0
		if candidates[0] != 0:
			candidates.append((self.p - candidates[0]))
		return candidates

	def from_y(self, y):
		# 1x^3 + 0x^2 + ax^1 + (b - y*y) = 0 mod p
		poly = [(self.b - y * y) % self.p, self.a, 0, 1]
		candidates = list(polyroots_prime(poly, self.p))
		assert len(candidates) != 0
		return candidates

	def check(self, pt):
		# does not check x val for origin(inf)
		rhs = ((pt.x**3) + (self.a * pt.x) + self.b - pt.y**2) % self.p
		return rhs == 0 or pt.z == 0

	def add(self, pt1, pt2):
		if pt1.is_origin():
			return pt2
		if pt2.is_origin():
			return pt1
		if (pt1.y + pt2.y) % self.p == 0 and pt1.x == pt2.x:
			return ec_point(0, 0)
		if pt1.x == pt2.x and pt1.y == pt2.y:
			temp = (((3 * pt1.x * pt1.x) + self.a) * modinv(2 * pt1.y, self.p)) % self.p
		else:
			temp = ((pt2.y - pt1.y) * modinv(pt2.x - pt1.x, self.p)) % self.p
		x = (temp * temp - pt1.x - pt2.x) % self.p
		y = (temp * (pt1.x - x) - pt1.y) % self.p

		pt = ec_point(x, y)
		print(pt)
		assert self.check(pt)
		return pt

	def invert(self, pt):
		return ec_point(pt.x, self.p - pt.y)

	def multiply(self, n, pt):
		if n == 0:
			return ec_point(0, 0)
		if n < 0:
			pt = self.invert(pt)
			assert self.check(pt)
			n = -n
		curr_bit_no = pt
		res = ec_point(0, 0)
		while n > 0:
			if n % 2 == 1:
				res = self.add(res, curr_bit_no)
			curr_bit_no = self.add(curr_bit_no, curr_bit_no)
			n = n // 2
		assert self.check(res)
		return res

	def mult(self, n, pt):
		return self.multiply(n, pt)

	def cardinality(self):
		if self.card is not None:
			return self.card
		if self.p < 200:
			self.card = self._cardinality_naive()
		else:
			self.card = self._cardinality_schoof()
		return self.card

	def _cardinality_naive(self):
		if self.card is not None:
			return self.card
		self.card = 0
		for i in range(self.p):
			self.card += len(self.from_x(i))
		return self.card

	def _cardinality_schoof(self):
		# TODO - needs math for endomorphism, Ltorsion, division polynomials
		raise NotImplementedError(
		    dedent("""
				Schoof algo is currently not implemented.
				In the meantime, use sage -
				sage_code(EllipticCurve(GF({self.p}), [{self.a},{self.b}]).order())
				OR
				add_sage_to_path();import sage.all as s; s.EllipticCurve(GF({self.p}), [{self.a},{self.b}]).order()"""))

	def _pollard_rho_step(self, r, a, b, g, f, n):
		# Hash fn = x mod 3
		choice = (r.x + r.y) % 3
		if choice == 0:
			return (self.multiply(2, r), (2 * a) % n, (2 * b) % n)
		elif choice == 1:
			return (self.add(f, r), a, b + 1)
		else:
			return (self.add(g, r), a + 1, b)

	def pollard_rho(self, g, f, n):
		# ri = ai*g + bi*f
		ri = g
		ai = 1
		bi = 0
		ri2 = g
		ai2 = 1
		bi2 = 0
		# One step for ri, 2 steps for ri2
		ri, ai, bi = self._pollard_rho_step(ri, ai, bi, g, f, n)
		ri2, ai2, bi2 = self._pollard_rho_step(ri2, ai2, bi2, g, f, n)
		ri2, ai2, bi2 = self._pollard_rho_step(ri2, ai2, bi2, g, f, n)

		while ri != ri2:
			ri, ai, bi = self._pollard_rho_step(ri, ai, bi, g, f, n)
			ri2, ai2, bi2 = self._pollard_rho_step(ri2, ai2, bi2, g, f, n)
			ri2, ai2, bi2 = self._pollard_rho_step(ri2, ai2, bi2, g, f, n)

		# Now, ri = ri2
		# ai*g + bi*f = ai2*g + bi2*f
		# (ai-ai2)*g = (bi2-bi)*f
		# x = (ai-ai2)/(bi2-bi) = (ai2-ai)/(bi-bi2)
		x = ((ai2 - ai) * modinv(bi - bi2, n)) % n
		return x
