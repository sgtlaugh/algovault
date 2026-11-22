"""
Discrete logarithm by baby-step giant-step

discrete_log(a, b, mod): smallest x >= 0 with pow(a, x, mod) == b % mod, None if no such x exists

Complexity: O(sqrt(mod)) modular multiplications and dict operations, O(sqrt(mod)) dict entries
"""

def discrete_log(a, b, mod):
    """
    returns smallest non-negative x s.t. pow(a, x, mod) == b % mod or None if no such x exists.
    Note: works even if a, b and mod are not pairwise coprime.
    a^0 = 1 (including 0^0), so x = 0 whenever b % mod == 1 % mod, which is every b when mod = 1
    """

    b %= mod
    if b == 1 % mod:
        return 0

    mp = {}
    e, n = 1, int(mod**0.5) + 2

    for i in range(n + 3):
        if e == b:
            return i

        mp[b * e % mod] = i
        e = e * a % mod

    v = e = pow(a, n, mod)
    for i in range(2, n + 3):
        e = e * v % mod
        if e in mp:
            x = n * i - mp[e]
            return x if pow(a, x, mod) == b else None


MOD = 100000007


def main():
    import sys

    data = sys.stdin.buffer.read().split()
    pos = 0
    t = int(data[pos])
    pos += 1
    out = []
    for case in range(1, t + 1):
        n, k, b, r = map(int, data[pos:pos + 4])
        pos += 4
        blocked = set()
        rows = 1
        for _ in range(b):
            x, y = int(data[pos]), int(data[pos + 1])
            pos += 2
            blocked.add((x, y))
            rows = max(rows, x)

        free = n - sum(1 for x, y in blocked if x == 1)
        free += sum(1 for x, y in blocked if x < rows and (x + 1, y) not in blocked)
        other = rows * n - b - free
        ways = pow(k, free, MOD) * pow(k - 1, other, MOD) % MOD

        if ways == r:
            out.append('Case %d: %d' % (case, rows))
            continue

        last = sum(1 for x, y in blocked if x == rows)
        ways = ways * pow(k, last, MOD) * pow(k - 1, n - last, MOD) % MOD
        step = pow(k - 1, n, MOD)
        res = discrete_log(step, r * pow(ways, MOD - 2, MOD), MOD)
        out.append('Case %d: %d' % (case, rows + 1 + res))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
