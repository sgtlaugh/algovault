"""
Lagrange Polynomial Interpolation
lagrange(x_values, y_values, nth_term): value at nth_term of the lowest degree polynomial through the points, exact

Complexity: O(n^2) big integer multiplications and n Fraction additions for n points, the integers grow with n
"""

from fractions import Fraction


def lagrange(x_values, y_values, nth_term):
    """
    If a polynomial has degree k, then it can be uniquely identified if it's values are known in k+1 distinct points

    :param x_values: X values of the polynomial
    :param y_values: Y values of the polynomial
    :param nth_term: nth_term of the polynomial
    :return: evaluates and returns the nth_term of the polynomial in O(k^2), given at least k+1 unique values were provided
             an int when the value is integral, a Fraction otherwise
    """

    if len(x_values) != len(y_values):
        raise ValueError('The X and Y values should be of the same length')

    res = Fraction(0)
    for i in range(len(x_values)):
        x, y = 1, 1
        for j in range(len(x_values)):
            if i != j:
                x *= (nth_term - x_values[j])
                y *= (x_values[i] - x_values[j])

        res += Fraction(x * y_values[i], y)  # a single term can be fractional even when the sum is not

    return res.numerator if res.denominator == 1 else res


def count(n):
    cells = [(r, c) for r in range(n) for c in range(n)]
    safe = [0] * len(cells)
    for i, (r, c) in enumerate(cells):
        for j, (r2, c2) in enumerate(cells):
            if j > i and r != r2 and c != c2 and r - c != r2 - c2 and r + c != r2 + c2:
                safe[i] |= 1 << j

    res = 0
    for i in range(len(cells)):
        for j in range(i + 1, len(cells)):
            if safe[i] >> j & 1:
                res += (safe[i] & safe[j]).bit_count()

    return res


def main():
    import sys

    n = int(sys.stdin.read())
    x_values = [n % 2 + 2 * i for i in range(1, 8)]
    y_values = [count(x) for x in x_values]
    print(lagrange(x_values, y_values, n))


if __name__ == '__main__':
    main()
