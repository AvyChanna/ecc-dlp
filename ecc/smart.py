from random import randrange

from ec import ec_curve
from utils import modinv


def hensel_lift(curve, P):
	x, y, _ = map(int, tuple(P))
	p = curve.p
	t = (((x * x * x + curve.a * x + curve.b) - y * y) // p) % p
	t = (t * modinv(2 * y, p)) % p
	return list(map(int, (x, y + (curve.p * t))))


def smart(curve, P, Q):
	A = curve.a
	x1, y1 = hensel_lift(curve, P)
	x2, y2 = hensel_lift(curve, Q)
	lifted_p = curve.p**2
	lifted_a = (y2**2 - y1**2 - (x2**3 - x1**3))
	lifted_a = (lifted_a * modinv(x2 - x1, lifted_p)) % lifted_p
	lifted_b = (y1**2 - x1**3 - A * x1) % lifted_p
	modulo = curve.p**2
	# do not verify curve params
	lifted_curve = ec_curve(lifted_a, lifted_b, lifted_p, verify=False)
	lifted_pt1 = lifted_curve.mult(curve.p - 1, lifted_curve(x1, y1))
	lifted_pt2 = lifted_curve.mult(curve.p - 1, lifted_curve(x2, y2))
	dx1 = ((lifted_pt1.x - x1) // curve.p) % modulo
	dy1 = ((lifted_pt2.x - x2) // curve.p) % modulo
	dx2 = lifted_pt1.y - y1
	dy2 = lifted_pt2.y - y2
	m = (dy1 * dx2 * modinv(dx1 * dy2, modulo)) % modulo
	return m % curve.p


def test_smart():
	a = 425706413842211054102700238164133538302169176474
	b = 203362936548826936673264444982866339953265530166
	p = 730750818665451459112596905638433048232067471723
	order = 730750818665451459112596905638433048232067471723
	k = randrange(0, order)
	curve = ec_curve(a, b, p)
	pt1 = curve(
	    282839918090522288605124127127354425085389489933,
	    575310508344796277762061084697250650618482060609,
	)
	pt2 = curve.multiply(k, pt1)
	res = smart(curve, pt1, pt2)
	assert res == k


if __name__ == "__main__":
	test_smart()
