from stress import rng, scaled

from berlekamp_massey import solve_linear_recurrence


def iterate(coefficients, start, count, mod):
    """f(i) = sum c[j] * f(i - 1 - j)"""
    terms = list(start)
    while len(terms) < count:
        terms.append(sum(c * terms[-1 - j] for j, c in enumerate(coefficients)) % mod)
    return terms


def main():
    for _ in range(scaled(150)):
        mod = rng.choice([3, 5, 7, 998244353, 10**9 + 7])
        k = rng.randint(1, 6)
        coefficients = [rng.randrange(mod) for _ in range(k)]
        terms = iterate(coefficients, [rng.randrange(mod) for _ in range(k)], 300, mod)

        base = terms[:2 * k + 2 * rng.randint(0, 3)]
        shifted = [x + mod * rng.randint(-3, 3) for x in base]  # the solver reduces inputs itself
        for _ in range(5):
            n = rng.randrange(300)
            assert solve_linear_recurrence(shifted, n, mod) == terms[n]

        # Too far to iterate: the terms around n must still satisfy the recurrence
        n = rng.randint(10**17, 10**18)
        window = [solve_linear_recurrence(base, n - j, mod) for j in range(k + 1)]
        assert window[0] == sum(c * window[1 + j] for j, c in enumerate(coefficients)) % mod


if __name__ == '__main__':
    main()
