from math import gcd, sqrt


def extended_gcd(a, b):
    if not b:
        return a, 1, 0

    g, y, x = extended_gcd(b, a % b)
    return g, x, y - ((a // b) * x)


def discrete_log(a, b, m):
    """returns smallest non-negative x s.t. pow(a, x, m) == b or None if no such x exists"""

    if b >= m:
        return None

    if b == 1 or m == 1:
        return 0

    i, t = 1, a
    while (1 << i) <= m:
        if t % m == b:
            return i
        i, t = i + 1, t * a % m

    raw_b, raw_m = b, m
    c, d, g = 0, 1, gcd(a, m)
    while g > 1:
        b, m = b // g, m // g
        c, d = c + 1, d * (a // g) % m
        g = gcd(a, m)

    pos = {}
    p, n = 1, int(sqrt(m - 0.666)) + 1
    for i in range(n):
        if p not in pos:
            pos[p] = i + 1
        p = p * a % m

    for i in range(n):
        g, x, y = extended_gcd(d, m)
        g = ((x * b // g) % m + m) % m
        if g in pos:
            x = i * n + pos[g] + c - 1
            if pow(a, x, raw_m) == raw_b:
                return x

        d = d * p % m

    return None


def main():
    assert discrete_log(1, 1, 1) is None
    assert discrete_log(2, 1, 3) == 0
    assert discrete_log(2, 3, 3) is None
    assert discrete_log(2, 0, 4) == 2
    assert discrete_log(6, 0, 8) == 3
    assert discrete_log(2, 3, 6) is None
    assert discrete_log(2, 6, 10) == 4
    assert discrete_log(5, 33, 58) == 9
    assert discrete_log(3589, 58, 97) is None
    assert discrete_log(3589, 1, 97) == 0
    assert discrete_log(1000000000, 666666667, 1000000007) == 942190576
    assert discrete_log(97383811, 35158600862, 100000000009) == 127906939


if __name__ == '__main__':
    main()
