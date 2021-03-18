from textwrap import wrap


def norm_hex(data: str) -> str:
	res = '0' * (len(data) % 2) + data
	if res == '00':
		return res
	while res.startswith('00'):
		res = res[2:]
	return res


def i2h(data: int):
	res = hex(int(data))[2:]
	return norm_hex(res)


def i2b(data: int):
	return h2b(i2h(data))


def i2a(data: int):
	return h2a(i2h(data))


def h2i(data: str):
	return int(data, 16)


def h2b(data: str):
	return bytes.fromhex(norm_hex(data))


def h2a(data: str):
	res = wrap(norm_hex(data), 2)
	return [int(i, 16) for i in res]


def b2i(data: bytes):
	return int.from_bytes(data, 'big')


def b2h(data: bytes):
	res = data.hex()
	return norm_hex(res)


def b2a(data: bytes):
	return h2a(b2h(data))


def a2i(data: list[int]):
	return h2i(a2h(data))


def a2h(data: list[int]):
	return "".join([i2h(i) for i in data])


def a2b(data: list[int]):
	while data[0] == 0:
		data = data[1:]
	return bytes(data)
