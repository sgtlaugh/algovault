from stress import rng, scaled

from discrete_log import discrete_log


def brute(a, b, mod):
    """The library's convention: with mod = 1, x = 0 matches only a literal b == 1 and every other b needs x = 1"""
    if b == 1:
        return 0
    b %= mod
    e = 1
    for x in range(mod + 1):
        if e == b:
            return x
        e = e * a % mod
    return None


def main():
    for it in range(scaled(5000)):
        mod = rng.randint(1, 60 if it % 10 else 3000)
        a = rng.randrange(mod) if it % 3 else rng.randint(-mod, 3 * mod)
        b = pow(a, rng.randint(0, 2 * mod), mod) if it % 2 else rng.randint(-3 * mod, 3 * mod)
        assert discrete_log(a, b, mod) == brute(a, b, mod), (a, b, mod)

    # Moduli too large for brute force: b = a^k, so a solution no larger than k must come back
    for _ in range(scaled(20)):
        mod = rng.randint(10**6, 10**10)
        a, k = rng.randrange(mod), rng.randrange(mod)
        b = pow(a, k, mod)
        x = discrete_log(a, b, mod)
        assert x is not None and 0 <= x <= k and pow(a, x, mod) == b


if __name__ == '__main__':
    main()
