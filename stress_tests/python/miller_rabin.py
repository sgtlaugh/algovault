from stress import rng, scaled

from miller_rabin import is_prime

SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41]


def reference(n):
    """Miller-Rabin with the first 13 primes as bases, deterministic below 3.3e24 and independent of the library's bases"""
    if n < 2:
        return False
    for p in SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in SMALL_PRIMES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def random_prime(bits):
    while True:
        n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if reference(n):
            return n


def main():
    limit = 300000
    sieve = [True] * limit
    sieve[0] = sieve[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = [False] * len(sieve[i * i::i])
    for n in range(-5, limit):
        assert is_prime(n) == (n >= 0 and sieve[n]), n

    # Strong pseudoprimes to many small bases and Carmichael numbers, the classic traps for weak base sets
    traps = [2047, 1373653, 25326001, 3215031751, 2152302898747, 3474749660383, 341550071728321, 3825123056546413051,
             318665857834031151167461, 561, 41041, 825265, 321197185, 5394826801, 232250619601, 9746347772161]
    for n in traps:
        assert not is_prime(n), n

    for _ in range(scaled(3000)):
        n = rng.randint(2, 2**64 - 1) if rng.randint(0, 1) else rng.randint(2, 10**12)
        assert is_prime(n) == reference(n), n

    # Beyond 2^64 the fixed bases are only probabilistic, but a prime is never rejected and a product of two large primes always is
    for _ in range(scaled(200)):
        p, q = random_prime(rng.randint(33, 90)), random_prime(rng.randint(33, 90))
        assert is_prime(p) and is_prime(q) and not is_prime(p * q) and not is_prime(p * p)


if __name__ == '__main__':
    main()
