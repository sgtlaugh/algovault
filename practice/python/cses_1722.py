"""
Fast Fibonacci

fibonacci(n, modulo=None): the n-th Fibonacci number by fast doubling, reduced mod modulo when given

Complexity:
  - with modulo: O(log n) multiplications of O(log modulo) bit integers
  - without: O(log n) multiplications, but F(n) has about 0.7 n bits, so the last ones dominate at O(n^1.585) with CPython's Karatsuba
"""

def fibonacci(n, modulo=None):
    """
    calculates the n-th fibonacci number using fast doubling algorithm
    fib(0) = 0, fib(1) = 1 and fib(n) = fib(n-1) + fib(n-2) for n > 1

    :param n: a non-negative integer
    :return: n-th fibonacci number
    """

    a, b = 0, 1
    for bit in bin(n)[2:-1]:
        a, b = a * (2 * b - a), a * a + b * b
        if bit == '1':
            a, b = b, a + b
        if modulo is not None:
            a, b = a % modulo, b % modulo

    # The last bit only needs F(n), skipping F(n + 1) saves the largest multiplication
    x = a * a + b * b if n & 1 else a * (2 * b - a)
    return x if modulo is None else x % modulo


def main():
    import sys

    n = int(sys.stdin.read())
    print(fibonacci(n, 10**9 + 7))


if __name__ == '__main__':
    main()
