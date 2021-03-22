from collections import namedtuple

from utils import isprime, modinv, polyroots_prime, sqrtmod_prime


class ec_point(namedtuple("ec_point", "x y z", defaults=(1, ))):
	def is_origin(self):
		return self.z == 0

	def __repr__(self):
		return f"pt({self.x}, {self.y})" if self.z == 1 else "Pt(Origin)"


class ec_curve:
	# Elliptic curve over finite prime field Z/nZ
	def __init__(self, a: int, b: int, p: int, order: int = None, verify: bool = True):
		self.p = p
		self.a = a % p
		self.b = b % p
		a4 = 4 * (self.a**3)
		b27 = 27 * (self.b**2)
		a4b27 = a4 + b27
		self.discriminant = -16 * a4b27
		if verify:
			assert isprime(p), f"{p} is not prime and verify=True in curve params"
			assert self.discriminant != 0, f"discriminant {self.discriminant} != 0 and verify=True in curve params"
		self.j_invariant = 1728 + (a4 * modinv(a4b27, self.p))
		self.card = order

	def __repr__(self):
		return f"elliptic curve y^2 = x^3 + {self.a}x + {self.b} over Zmod({self.p})"

	def from_x(self, x):
		# y2 = [x^3 + ax +b] mod p
		y2 = (x * x * x + self.a * x + self.b) % self.p
		# screw tonelli shanks, hensel lift rocks
		candidate = sqrtmod_prime(y2, self.p)
		if (candidate**2) % self.p == y2:
			return list(set([ec_point(x, candidate), ec_point(x, (self.p - candidate) % self.p)]))
		# for convenience, this won't raise
		return []

	def from_y(self, y):
		# Find solution(s) to x^3 + ax + (b - y^2) = 0 mod p
		poly = [(self.b - y * y) % self.p, self.a, 0, 1]
		candidates = list(polyroots_prime(poly, self.p))
		assert len(candidates) != 0, "no point with y={y} on {self}"
		return candidates

	def check(self, pt):
		# does not check xy val for origin(inf)
		rhs = ((pt.x**3) + (self.a * pt.x) + self.b - (pt.y**2)) % self.p
		return rhs == 0 or pt.is_origin()

	def assert_check(self, pt):
		assert self.check(pt), f"{pt} does not lie on {self}"

	def add(self, pt1, pt2):
		if pt1.is_origin():
			return pt2
		if pt2.is_origin():
			return pt1
		if (pt1.y + pt2.y) % self.p == 0 and pt1.x == pt2.x:
			return self.origin()
		if pt1.x == pt2.x and pt1.y == pt2.y:
			temp = (((3 * pt1.x * pt1.x) + self.a) * modinv(2 * pt1.y, self.p)) % self.p
		else:
			temp = ((pt2.y - pt1.y) * modinv(pt2.x - pt1.x, self.p)) % self.p
		x = (temp * temp - pt1.x - pt2.x) % self.p
		y = (temp * (pt1.x - x) - pt1.y) % self.p

		pt = ec_point(x, y)
		self.assert_check(pt)
		return pt

	def invert(self, pt):
		return ec_point(pt.x, self.p - pt.y)

	def multiply(self, n, pt):
		if n == 0:
			return self.origin()
		if n < 0:
			pt = self.invert(pt)
			self.assert_check(pt)
			n = -n
		curr_bit_no = pt
		res = self.origin()
		while n > 0:
			if n % 2 == 1:
				res = self.add(res, curr_bit_no)
			curr_bit_no = self.add(curr_bit_no, curr_bit_no)
			n = n // 2
		self.assert_check(res)
		return res

	def mult(self, n, pt):
		return self.multiply(n, pt)

	def mul(self, n, pt):
		return self.multiply(n, pt)

	def cardinality(self):
		if self.card is not None:
			return self.card
		if self.p < 5000:
			self.card = self._cardinality_naive()
		else:
			self.card = self._cardinality_schoof()
		return self.card

	def _cardinality_naive(self):
		if self.card is not None:
			return self.card
		self.card = 0
		for i in range(self.p):
			pts = self.from_x(i)
			if len(pts) != 0:
				self.card += len(pts)
		return self.card

	def _cardinality_schoof(self):
		# TODO - needs math for endomorphism, torsion, division polynomials
		raise NotImplementedError("Schoof algorithm is currently not implemented. In the meantime, use sage")

	def __call__(self, x, y, z=1):
		res = ec_point(x, y, z)
		self.assert_check(res)
		return res

	@staticmethod
	def origin():
		return ec_point(0, 1, 0)
