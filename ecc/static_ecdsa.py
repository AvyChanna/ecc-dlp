from collections import namedtuple
from hashlib import new as new_hash
from typing import Callable

from ec import ec_curve, ec_point
from utils import b2i, i2b, modinv

Public_key = namedtuple("Public_key", "gen order curve_pt")
Private_key = namedtuple("Private_key", "gen order secret")
Signature = namedtuple("Signature", "r s")


def ecdsa_gen_keys(curve: ec_curve, gen: ec_point, order: int, secret: int):
	curve_pt = curve.mult(secret, gen)
	public_key = Public_key(gen, order, curve_pt)
	private_key = Private_key(gen, order, secret)
	return public_key, private_key


def ecdsa_verify(curve: ec_curve, pub_key: Public_key, msg_hash: bytes, signature: tuple):
	gen, order, curve_pt = pub_key
	r, s = signature

	if r < 1 or r > order - 1:
		return False
	if s < 1 or s > order - 1:
		return False

	c = modinv(s, order)
	a1 = (msg_hash * c) % order
	a2 = (r * c) % order
	xy = curve.add(curve.mult(a1, gen), curve.mult(a2, curve_pt))
	ver = xy.x % order
	return r == ver


def ecdsa_sign(curve: ec_curve, priv_key: Private_key, msg_hash: int, nonce: int):
	gen, order, secret = priv_key
	p = curve.mult(nonce, gen)
	r = p.x % order
	if r == 0:
		raise ValueError("Got r=0. Unlucky.")
	s = (modinv(nonce, order) * (msg_hash + (secret * r) % order)) % order
	if s == 0:
		raise ValueError("Got s=0. Unlucky.")
	return (r, s)


def ecdsa_sign_static(curve: ec_curve, priv_key: Private_key, msg_hash: int, hash_func: Callable[[int, int], int]):
	nonce = hash_func(priv_key.secret, msg_hash)  # nonce
	return ecdsa_sign(curve, priv_key, msg_hash, nonce)


def insecure_hash(secret: int, msg_hash: int, hash_name: str = "md5") -> int:
	# Here insecure means that the hash function is partially/fully insecure, due to
	# - weak/strong collisions
	# - length extension attacks
	# - fastcoll/unicoll type attacks
	# and so on.
	secret_bytes = i2b(secret)
	msg_hash_bytes = i2b(msg_hash)
	hash_func = lambda x: new_hash(hash_name, x).digest()
	return b2i(hash_func(secret_bytes + msg_hash_bytes))


def hmac(key: bytes, msg: bytes, hash_func: Callable, block_size: int):
	if len(key) > block_size:
		key = hash_func(key)

	key = key + b'\x00' * (block_size - len(key))

	o_key_pad = bytes(a ^ b for a, b in zip(key, b'\x5c' * block_size))
	i_key_pad = bytes(a ^ b for a, b in zip(key, b'\x36' * block_size))

	return hash_func(o_key_pad + hash_func(i_key_pad + msg))


def secure_hash(secret: int, msg_hash: int, hash_name: str = "md5") -> int:
	secret_bytes = i2b(secret)
	msg_hash_bytes = i2b(msg_hash)
	hash_func = lambda x: new_hash(hash_name, x).digest()
	hash_block_size = new_hash(hash_name).block_size
	return b2i(hmac(secret_bytes, msg_hash_bytes, hash_func, hash_block_size))


def get_test_curve256():
	# This is a "secure enough" curve.
	p = 115792089210356248762697446949407573530086143415290314195533631308867097853951
	a = 115792089210356248762697446949407573530086143415290314195533631308867097853948  # -3 mod p
	b = 41058363725152142129326129780047268409114441015993725554835256314039467401291
	order = 115792089210356248762697446949407573529996955224135760342422259061068512044369
	curve = ec_curve(a, b, p, order)

	Gx = 48439561293906451759052585252797914202762949526041747995844080717082404635286
	Gy = 36134250956749795798585127919587881956611106672985015071877198253568414405109
	G = curve(Gx, Gy)
	return curve, G, order
