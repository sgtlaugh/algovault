def is_prime(n):
    _BASES = [2, 450775, 1795265022, 9780504, 28178, 9375, 325]

    if n < 5 or n & 1 == 0 or n % 3 == 0:
        return n == 2 or n == 3

    s = ((n - 1) & (1 - n)).bit_length() - 1
    return not any(_is_composite(a, n, s) for a in _BASES)


def _is_composite(a, n, s):
    p = pow(a, n >> s, n)
    if p == 1 or p == n - 1 or a % n == 0:
        return False

    for _ in range(s):
        p = p * p % n
        if p == n - 1:
            return False

    return True


def main():
    assert is_prime(1) is False
    assert is_prime(2) is True
    assert is_prime(3) is True
    assert is_prime(4) is False
    assert is_prime(97) is True
    assert is_prime(299210837) is True
    assert is_prime(1000000007) is True
    assert is_prime(2147483647) is True
    assert is_prime(1000000007 * 2147483647) is False
    assert is_prime(1000003 * 474119 * 699053) is False
    assert is_prime(4547337172376300111955330758342147474062293202868155909489) is True


if __name__ == '__main__':
    main()
