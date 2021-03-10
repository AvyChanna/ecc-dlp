from random import randrange

from bsgs import bsgs
from ec import ec_curve
from utils import crt, factorint


def pohlig(curve, g, f, n, limit=0):
	subgroups = n
	if isinstance(n, int):
		subgroups = factorint(n)
	factors = [i**j for i, j in subgroups.items()]
	prod = 1
	for i in factors:
		prod *= i
	exponents = [prod // i for i in factors]
	remainders = []
	for factor, power in zip(factors, exponents):
		if limit != 0 and factor > limit:
			continue
		g_power_k = curve.mult(power, g)
		assert curve.check(g_power_k)
		r_power_k = curve.mult(power, f)
		remainders.append(bsgs(curve, g_power_k, r_power_k, factor))
	return crt(remainders, factors)


def test_pohlig():
	a = 2
	b = 3
	p = 1125899839733759
	order = 1125899867612160
	k = randrange(2, order)
	curve = ec_curve(a, b, p)
	pt1 = curve(436757568245484, 726713018309225)
	pt2 = curve.multiply(k, pt1)
	res = pohlig(curve, pt1, pt2, order)
	assert res == k


if __name__ == "__main__":
	test_pohlig()
