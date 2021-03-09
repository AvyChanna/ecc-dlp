import math


def bsgs(curve, g, f, n):
	m = math.ceil(math.sqrt(n))
	baby_step = {curve.multiply(i, g): i for i in range(m)}
	m_inv = curve.multiply(-m, g)
	y = f
	for i in range(m):
		if y in baby_step:
			return i * m + baby_step[y]
		y = curve.add(f, m_inv)
	return None
