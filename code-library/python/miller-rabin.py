from random import randint


def miller_rabin(n, trials=20):
    def _is_composite(a, d, n, s):
        if pow(a, d, n) == 1:
            return False

        for i in range(s):
            if pow(a, 2 ** i * d, n) == n - 1:
                return False

        return True

    if n < 2:
        return False

    small_primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}
    if any(n % p == 0 for p in small_primes):
        return n in small_primes

    x = str(bin(n - 1))
    s = len(x) - len(x.rstrip('0'))
    d = (n - 1) >> s

    witness_list = {randint(2, n - 1) for _ in range(trials)}
    return not any(_is_composite(a, d, n, s) for a in witness_list)


def main():
    print(miller_rabin(2))                                                           # True
    print(miller_rabin(3))                                                           # True
    print(miller_rabin(4))                                                           # False
    print(miller_rabin(97))                                                          # True
    print(miller_rabin(1000000007))                                                  # True
    print(miller_rabin(2147483647))                                                  # True
    print(miller_rabin(1000000007 * 2147483647))                                     # False
    print(miller_rabin(4547337172376300111955330758342147474062293202868155909489))  # True


if __name__ == '__main__':
    main()
