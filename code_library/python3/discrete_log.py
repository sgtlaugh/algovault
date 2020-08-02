def discrete_log(a, b, mod):
    """
    returns smallest non-negative x s.t. pow(a, x, mod) == b % mod
    or None if no such x exists
    """

    mp = {}
    e, n = 1, int(mod**0.5) + 2

    for j in range(n + 3):
        if e == b:
            return j

        mp[b * e % mod] = j
        e = e * a % mod

    v = e = pow(a, n, mod)
    for i in range(2, n + 3):
        e = e * v % mod
        if e in mp:
            x = n * i - mp[e]
            return x if pow(a, x, mod) == b else None


def main():
    assert discrete_log(0, 0, 1) == 1
    assert discrete_log(1, 1, 1) == 0
    assert discrete_log(2, 1, 3) == 0
    assert discrete_log(2, 3, 3) is None
    assert discrete_log(2, 0, 4) == 2
    assert discrete_log(6, 0, 8) == 3
    assert discrete_log(6, 8, 16) == 3
    assert discrete_log(2, 3, 6) is None
    assert discrete_log(2, 6, 10) == 4
    assert discrete_log(5, 33, 58) == 9
    assert discrete_log(3589, 58, 97) is None
    assert discrete_log(3589, 1, 97) == 0
    assert discrete_log(1000000000, 666666667, 1000000007) == 942190576
    assert discrete_log(97383811, 35158600862, 100000000009) == 127906939


if __name__ == '__main__':
    main()
