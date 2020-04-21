from random import randint


def miller_rabin(n, ntrials=15):
    _SMALL_PRIMES = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}

    if n < 2:
        return False

    if any(n % p == 0 for p in _SMALL_PRIMES):
        return n in _SMALL_PRIMES

    s = str(bin(n - 1))
    s = len(s) - len(s.rstrip('0'))

    witness_list = {randint(2, n - 1) for _ in range(ntrials)}
    return not any(_is_composite(a, n, s) for a in witness_list)


def _is_composite(a, n, s):
    d = (n - 1) >> s
    if pow(a, d, n) == 1:
        return False

    return all(pow(a, (1 << i) * d, n) != (n - 1) for i in range(s))


def main():
    print(miller_rabin(2))  # True
    print(miller_rabin(3))  # True
    print(miller_rabin(4))  # False
    print(miller_rabin(97))  # True
    print(miller_rabin(1000000007))  # True
    print(miller_rabin(2147483647))  # True
    print(miller_rabin(1000000007 * 2147483647))  # False
    print(miller_rabin(4547337172376300111955330758342147474062293202868155909489))  # True


if __name__ == '__main__':
    main()
