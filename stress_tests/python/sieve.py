from stress import rng, scaled

from sieve import sieve, primes_up_to


def is_prime(x):
    """Trial division, independent of the sieve"""
    if x < 2:
        return False
    d = 2
    while d * d <= x:
        if x % d == 0:
            return False
        d += 1
    return True


def main():
    for n in range(0, 300):
        flags = sieve(n)
        assert len(flags) == n + 1
        assert all(flags[x] == is_prime(x) for x in range(n + 1)), n

    for _ in range(scaled(20)):
        n = rng.randint(1000, 200000)
        flags = sieve(n)
        for _ in range(500):
            x = rng.randint(0, n)
            assert flags[x] == is_prime(x), (n, x)
        assert primes_up_to(n)[-1] == max(x for x in range(n - 200, n + 1) if is_prime(x))

    assert len(primes_up_to(10000000)) == 664579


if __name__ == '__main__':
    main()
