import math
from random import randrange

from ec import ec_curve


def bsgs(curve, g, f, n):
	m = math.ceil(math.sqrt(n))
	baby_step = {curve.multiply(i, g): i for i in range(m)}
	m_inv = curve.multiply(-m, g)
	y = f
	for i in range(m):
		if y in baby_step:
			return i * m + baby_step[y]
		y = curve.add(y, m_inv)
	return None


def test_bsgs():
	a = 2
	b = 3
	p = 4111
	order = 4120
	k = randrange(2, order)
	curve = ec_curve(a, b, p)
	pt1 = curve(2672, 2565)
	pt2 = curve.mult(k, pt1)
	res = bsgs(curve, pt1, pt2, order)
	assert res == k


if __name__ == "__main__":
	test_bsgs()
