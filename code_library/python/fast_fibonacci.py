import sys


def fibonacci(n, modulo=None):
    """
    calculates the n-th fibonacci number using fast doubling algorithm
    fib(0) = 0, fib(1) = 1 and fib(n) = fib(n-1) + fib(n-2) for n > 1

    :param n: a non-negative integer
    :return: n-th fibonacci number
    """

    def _fib(n, modulo):
        if n == 0:
            return 0, 1

        a, b = _fib(n >> 1, modulo)
        x, y = a * (2 * b - a), a * a + b * b

        if n & 1:
            x, y = y, x + y

        if modulo is not None:
            x, y = x % modulo, y % modulo

        return x, y

    x, y = _fib(n, modulo)
    return x if modulo is None else x % modulo


def main():
    sys.setrecursionlimit(2 ** 22)
    
    assert fibonacci(0) == 0
    assert fibonacci(1) == 1
    assert fibonacci(2) == 1
    assert fibonacci(10) == 55
    assert fibonacci(100) == 354224848179261915075
    assert fibonacci(400) == 176023680645013966468226945392411250770384383304492191886725992896575345044216019675

    modulo = 10 ** 9 + 7
    assert fibonacci(0, modulo) == 0
    assert fibonacci(1, modulo) == 1
    assert fibonacci(2, modulo) == 1
    assert fibonacci(10, modulo) == 55
    assert fibonacci(100, modulo) == 687995182
    assert fibonacci(400, modulo) == 967250938
    assert fibonacci(10**1000, modulo) == 552179166
    
    assert len(str(fibonacci(10000))) == 2090
    assert len(str(fibonacci(100000))) == 20899
    assert len(str(fibonacci(1000000))) == 208988


if __name__ == '__main__':
    main()
