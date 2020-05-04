import sys


def fibonacci(n, modulo=None):
    """
    calculates the nth fibonacci number using fast doubling algorithm
    fib(0) = 0, fib(1) = 1 and fib(n) = fib(n-1) + fib(n-2) for n > 1

    :param n:
    :return: nth fibonacci number
    """

    def _fib(n):
        if n == 0:
            return 0, 1

        a, b = _fib(n >> 1)
        x, y = a * (2 * b - a), a * a + b * b

        if n & 1:
            x, y = y, x + y

        return x, y

    def _fibmod(n, modulo):
        if n == 0:
            return 0 % modulo, 1 % modulo

        a, b = _fibmod(n >> 1, modulo)
        x, y = a * (2 * b - a), a * a + b * b

        if n & 1:
            x, y = y, x + y

        return x % modulo, y % modulo

    x, y = _fib(n) if modulo is None else _fibmod(n, modulo)
    return x


def main():
    sys.setrecursionlimit(2 ** 22)

    print(fibonacci(0))  # 0
    print(fibonacci(1))  # 1
    print(fibonacci(2))  # 1
    print(fibonacci(10))  # 55
    print(fibonacci(100))  # 354224848179261915075
    print(fibonacci(400))  # 176023680645013966468226945392411250770384383304492191886725992896575345044216019675

    modulo = 10 ** 9 + 7
    print(fibonacci(0, modulo))  # 0
    print(fibonacci(1, modulo))  # 1
    print(fibonacci(2, modulo))  # 1
    print(fibonacci(10, modulo))  # 55
    print(fibonacci(100, modulo))  # 687995182
    print(fibonacci(400, modulo))  # 967250938
    print(fibonacci(10 ** 18, modulo))  # 209783453
    print(fibonacci(10 ** 1000, modulo))  # 552179166

    # the following works pretty fast too, ~1s locally
    # print(fibonacci(10000))
    # print(fibonacci(100000))
    # print(fibonacci(1000000))


if __name__ == '__main__':
    main()
