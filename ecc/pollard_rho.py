from math import gcd
from random import randrange

from ec import ec_curve
from utils import xgcd


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


def pollard_rho(curve, g, f, n):
	# ri1 = ai1*g + bi1*f
	ri1, ai1, bi1 = g, 1, 0
	ri2, ai2, bi2 = g, 1, 0

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
	import math
	print((ai2 - ai1))
	print((bi1 - bi2))
	lhs = (ai2 - ai1) % n
	rhs = (bi1 - bi2) % n
	print(lhs, rhs)
	x = (lhs * xgcd(rhs, n)[1]) % n
	# TODO
	return x


def test_pollard_rho():
	a = 45181635
	b = 124806060
	p = 1035418103
	order = 1035393330
	k = 884158742
	# k = randrange(2, order)
	curve = ec_curve(a, b, p)
	pt1 = curve(734830775, 908510221)
	pt2 = curve.mult(k, pt1)
	res = pollard_rho(curve, pt1, pt2, order)
	print(k, res)


test_pollard_rho()

# lhs = 867599962
# rhs = 513486386
# k = 884158742
# order = 1035393330
#
#
#
#
#
