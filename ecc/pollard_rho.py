from random import randrange

from ec import ec_curve
from utils import modinv


def pollard_rho_step(curve, g, f, n, r, a, b):
	# Hash fn = x mod 3, any hash that divides domain
	# into 3 nearly equal parts is acceptable
	choice = (r.x + r.y) % 3
	if choice == 0:
		return (curve.multiply(2, r), (2 * a) % n, (2 * b) % n)
	elif choice == 1:
		return (curve.add(f, r), a, b + 1)
	else:
		return (curve.add(g, r), a + 1, b)


def pollard_rho(curve, g, f, n, init_a=1, init_b=0, max_tries=3):
	if max_tries == 0:
		return None
	# ri1 = ai1*g + bi1*f
	ag = curve.mult(init_a, g)
	bf = curve.mult(init_b, f)
	ag_bf = curve.add(ag, bf)
	ri1, ai1, bi1 = ag_bf, init_a, init_b
	ri2, ai2, bi2 = ag_bf, init_a, init_b

	# One step for ri1, 2 steps for ri2
	ri1, ai1, bi1 = pollard_rho_step(curve, g, f, n, ri1, ai1, bi1)
	ri2, ai2, bi2 = pollard_rho_step(curve, g, f, n, ri2, ai2, bi2)
	ri2, ai2, bi2 = pollard_rho_step(curve, g, f, n, ri2, ai2, bi2)

	while ri1 != ri2:
		ri1, ai1, bi1 = pollard_rho_step(curve, g, f, n, ri1, ai1, bi1)
		ri2, ai2, bi2 = pollard_rho_step(curve, g, f, n, ri2, ai2, bi2)
		ri2, ai2, bi2 = pollard_rho_step(curve, g, f, n, ri2, ai2, bi2)

	# Now, ri1 = ri2
	# ai1*g + bi1*f = ai2*g + bi2*f
	# (ai1-ai2)*g = (bi2-bi1)*f
	# x = (ai1-ai2)/(bi2-bi1) = (ai2-ai1)/(bi1-bi2)

	# An edge case where we get 0*G = 0*F. In this case we start from scratch,
	# but with random(different) initial values. The probability that this will
	# continuously fail should be fairly low. Nonetheless, we exit after `max_tries` failures
	if (bi1 - bi2) % n == 0:
		return pollard_rho(curve, g, f, n, init_a + randrange(0, n), init_b + randrange(0, n), max_tries - 1)

	x = ((ai2 - ai1) * modinv(bi1 - bi2, n)) % n
	return x


def test_pollard_rho():
	a = 2
	b = 9
	p = 1035418103
	order = 1035356653
	k = randrange(2, order)
	curve = ec_curve(a, b, p)
	pt1 = curve(769278016, 752868328)
	pt2 = curve.mult(k, pt1)
	res = pollard_rho(curve, pt1, pt2, order)
	assert res == k


test_pollard_rho()
