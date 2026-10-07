from fractions import Fraction

from stress import rng, scaled

from lagrange_polynomial_interpolation import lagrange


def evaluate(coefficients, x):
    return sum(c * x**i for i, c in enumerate(coefficients))


def main():
    for it in range(scaled(2000)):
        degree = rng.randint(0, 8 if it % 10 else 40)
        coefficients = [rng.randint(-10**6, 10**6) for _ in range(degree + 1)]
        if it % 3 == 0:
            coefficients = [Fraction(c, rng.randint(1, 50)) for c in coefficients]  # values need not be integers

        x_values = rng.sample(range(-200, 200), degree + 1 + rng.randint(0, 3))  # distinct, unordered, gaps
        y_values = [evaluate(coefficients, x) for x in x_values]

        for _ in range(5):
            n = rng.choice([rng.randint(-300, 300), rng.randint(-10**18, 10**18), rng.choice(x_values)])
            res, expected = lagrange(x_values, y_values, n), evaluate(coefficients, n)
            assert res == expected and isinstance(res, int) == (Fraction(expected).denominator == 1)


if __name__ == '__main__':
    main()
