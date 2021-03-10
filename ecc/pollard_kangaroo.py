from math import ceil, sqrt

from ec import ec_curve


class PRF:
	def __init__(self, a, b):
		self.state = ceil(sqrt((b - a)) / 2)

	def get(self, P):
		if P.is_origin():
			return 1
		return int(P.x) % self.state + int(P.y) % self.state + 1


def pollard_kangaroo(G, curve, F, a, b, N=None):
	if N is None:
		N = ceil(sqrt(b - a))
	hash_func = PRF(a, b)

	# tame search
	x_tame = 0
	y_tame = curve.mult(b, G)

	while N > 0:
		x_tame += hash_func.get(y_tame)
		y_tame = curve.add(y_tame, curve.mult(hash_func.get(y_tame), G))
		N -= 1

	assert y_tame == curve.mult(b + x_tame, G)

	# wild search
	x_wild = 0
	y_wild = F

	upper_limit = b - a + x_tame
	while x_wild < upper_limit:
		x_wild += hash_func.get(y_wild)
		y_wild = curve.add(y_wild, curve.mult(hash_func.get(y_wild), G))

		if y_wild == y_tame:
			return b + x_tame - x_wild

	# index not found
	return None


def test_pollard_kangaroo():
	a = 0
	b = 7
	p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
	order = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
	k = 1000
	curve = ec_curve(a, b, p, verify=False)
	G = curve(
	    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
	    0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8,
	)
	F = curve.mult(k, G)
	res = pollard_kangaroo(G, curve, F, 0, 2000)
	assert res == k


if __name__ == "__main__":
	test_pollard_kangaroo()
