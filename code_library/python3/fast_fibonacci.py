def fibonacci(n):
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

    x, y = _fib(n)
    return x


def main():
    print(fibonacci(0))     # 0
    print(fibonacci(1))     # 1
    print(fibonacci(2))     # 1
    print(fibonacci(10))    # 55
    print(fibonacci(100))   # 354224848179261915075
    print(fibonacci(400))   # 176023680645013966468226945392411250770384383304492191886725992896575345044216019675

    # the following works pretty fast too, ~1s locally
    # print(fibonacci(10000))
    # print(fibonacci(100000))
    # print(fibonacci(1000000))


if __name__ == '__main__':
    main()
