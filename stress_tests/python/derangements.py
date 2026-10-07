from fractions import Fraction
from itertools import permutations
from math import comb, exp, factorial, perm

from stress import rng, scaled

from derangements import derangement, derangement_probability


def inclusion_exclusion(n, m):
    """Seatings of n people in m seats avoiding the k chosen fixed points, summed with alternating signs"""
    return sum((-1) ** k * comb(n, k) * perm(m - k, n - k) for k in range(n + 1)) if n <= m else 0


def main():
    for n in range(7):
        for m in range(8):
            brute = sum(all(seat != person for person, seat in enumerate(chosen)) for chosen in permutations(range(m), n))
            assert derangement(n, m) == brute, (n, m)

    for _ in range(scaled(500)):
        n = rng.randint(0, 120)
        m = n + rng.randint(-3, 60)
        assert derangement(n, m) == inclusion_exclusion(n, m)

    for n in range(60):
        exact = Fraction(derangement(n, n), factorial(n)) if n <= 20 else exp(-1)
        assert abs(derangement_probability(n) - float(exact)) < 1e-15
        assert n < 21 or abs(float(Fraction(inclusion_exclusion(n, n), factorial(n))) - exp(-1)) < 1e-15  # past 20 the true ratio matches 1 / e to double precision


if __name__ == '__main__':
    main()
