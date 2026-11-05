"""
Miller-Rabin Primality Test

is_prime(n): Miller-Rabin test of n against 7 fixed bases
The bases are proven deterministic only for n < 2^64, above that a True is probabilistic (a composite can pass)

Complexity: O(log n) modular multiplications of log n bit integers per base, 7 bases
"""

def is_prime(n):
    _BASES = [2, 450775, 1795265022, 9780504, 28178, 9375, 325]

    if n < 5 or n & 1 == 0 or n % 3 == 0:
        return n == 2 or n == 3

    s = ((n - 1) & (1 - n)).bit_length() - 1
    return not any(_is_composite(a, n, s) for a in _BASES)


def _is_composite(a, n, s):
    p = pow(a, n >> s, n)
    if p == 1 or p == n - 1 or a % n == 0:
        return False

    for _ in range(s - 1):
        p = p * p % n
        if p == n - 1:
            return False

    return True


def main():
    import sys

    data = sys.stdin.read().split()
    out = []
    for token in data[1:int(data[0]) + 1]:
        n = int(token) + 1
        while not is_prime(n):
            n += 1
        out.append(str(n))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
