from utils import crt, factorint


def pohlig(curve, g, f, n):
	subgroups = n
	if isinstance(n, int):
		subgroups = factorint(n)
	factors = [i * j for i, j in subgroups.items()]
	prod = 1
	for i in factors:
		prod *= i
	exponents = [prod // i for i in factors]
	remainders = []
	for factor, power in zip(factors, exponents):
		g_power_k = curve.multiply(power, g)
		r_power_k = curve.multiply(power, f)
		remainders.append(curve.bsgs(g_power_k, r_power_k, factor))
	return crt(remainders, factors)
