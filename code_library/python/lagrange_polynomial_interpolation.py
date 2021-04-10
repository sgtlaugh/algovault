def lagrange(x_values, y_values, nth_term):
    """
    If a polynomial has degree k, then it can be uniquely identified if it's values are known in k+1 distinct points

    :param x_values: X values of the polynomial
    :param y_values: Y values of the polynomial
    :param nth_term: nth_term of the polynomial
    :return: evaluates and returns the nth_term of the polynomial in O(k^2), given at least k+1 unique values were provided
    """

    if len(x_values) != len(y_values):
        raise ValueError('The X and Y values should be of the same length')

    res = 0
    for i in range(len(x_values)):
        x, y = 1, 1
        for j in range(len(x_values)):
            if i != j:
                x *= (nth_term - x_values[j])
                y *= (x_values[i] - x_values[j])

        res += (x * y_values[i] // y)

    return res


def main():
    """

    Example polyonmial - https://mathworld.wolfram.com/CircleDivisionbyChords.html
    f(x) = binomial(n, 4) + binomial(n, 2) + 1
    or, f(x) = (n^4 - 6n^3 + 23n^2 -18n + 24) / 24

    """

    x_values = [0, 1, 2, 3, 4]
    y_values = [1, 1, 2, 4, 8]

    assert lagrange(x_values, y_values, 0) == 1
    assert lagrange(x_values, y_values, 1) == 1
    assert lagrange(x_values, y_values, 4) == 8
    assert lagrange(x_values, y_values, 7) == 57
    assert lagrange(x_values, y_values, 10) == 256
    assert lagrange(x_values, y_values, 20) == 5036
    assert lagrange(x_values, y_values, 1000000000) == 41666666416666667624999999250000001


if __name__ == '__main__':
    main()
