from random import randint


def miller_rabin(n, ntrials=15):
    _SMALL_PRIMES = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}

    if n < 2:
        return False

    if any(n % p == 0 for p in _SMALL_PRIMES):
        return n in _SMALL_PRIMES

    s = ((n - 1) & (1 - n)).bit_length() - 1
    witness_list = {randint(2, n - 1) for _ in range(ntrials)}
    return not any(_is_composite(a, n, s) for a in witness_list)


def _is_composite(a, n, s):
    d = (n - 1) >> s
    if pow(a, d, n) == 1:
        return False

    return all(pow(a, (1 << i) * d, n) != (n - 1) for i in range(s))


def main():
    assert miller_rabin(2) is True
    assert miller_rabin(3) is True
    assert miller_rabin(4) is False
    assert miller_rabin(97) is True
    assert miller_rabin(1000000007) is True
    assert miller_rabin(2147483647) is True
    assert miller_rabin(1000000007 * 2147483647) is False
    assert miller_rabin(4547337172376300111955330758342147474062293202868155909489) is True


if __name__ == '__main__':
    main()
