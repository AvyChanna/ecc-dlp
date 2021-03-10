from bsgs import test_bsgs
from pohlig import test_pohlig
from pollard_kangaroo import test_pollard_kangaroo
from pollard_rho import test_pollard_rho
from smart import test_smart

tests = {
    "bsgs": test_bsgs,
    "pohlig": test_pohlig,
    "smart": test_smart,
    "pollard_rho": test_pollard_rho,
    "pollard_kangaroo": test_pollard_kangaroo,
}


def test():
	for name, test in tests.items():
		try:
			test()
			print(f"[*] PASSED TEST {name}")
		except:
			print(f"[!] FAILED TEST {name}")


if __name__ == "__main__":
	test()
