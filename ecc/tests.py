from bsgs import test_bsgs
from pohlig import test_pohlig
from pollard_kangaroo import test_pollard_kangaroo
from pollard_rho import test_pollard_rho
from smart import test_smart

tests = {
    "BABY_STEP_GIANT_STEP": test_bsgs,
    "POHLIG_HELLMAN": test_pohlig,
    "SMART": test_smart,
    "POLLARD_RHO": test_pollard_rho,
    "POLLARD_KANGAROO": test_pollard_kangaroo,
}


def test_all():
	for name, test in tests.items():
		try:
			test()
			print(f"\033[92m[*] PASSED TEST \033[94m{name}\033[0m")
		except KeyboardInterrupt:
			break
		except:  # pylint: disable=bare-except
			print(f"\033[91m[!] FAILED TEST \033[93m{name}")


if __name__ == "__main__":
	test_all()
