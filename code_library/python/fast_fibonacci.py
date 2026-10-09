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
    
    # Digit counts checked by exact comparison, str() refuses ints over 4300 digits since Python 3.11
    has_digits = lambda x, d: 10 ** (d - 1) <= x < 10 ** d
    assert has_digits(fibonacci(10000), 2090)
    assert has_digits(fibonacci(100000), 20899)
    assert has_digits(fibonacci(1000000), 208988)


if __name__ == '__main__':
    main()
