from stress import rng, scaled

from fast_fibonacci import fibonacci


def main():
    fib = [0, 1]
    while len(fib) < 3000:
        fib.append(fib[-1] + fib[-2])

    for n in range(3000):
        assert fibonacci(n) == fib[n]

    for _ in range(scaled(2000)):
        n, modulo = rng.randrange(3000), rng.choice([1, 2, 10, 97, 10**9 + 7, 2**61 - 1, rng.randint(1, 10**30)])
        assert fibonacci(n, modulo) == fib[n] % modulo

    # Too large to iterate: F(a + b) = F(a) F(b + 1) + F(a - 1) F(b)
    for _ in range(scaled(500)):
        modulo = rng.choice([10**9 + 7, rng.randint(2, 10**18)])
        a, b = rng.randint(1, 10**100), rng.randint(0, 10**100)
        f = lambda n: fibonacci(n, modulo)
        assert f(a + b) == (f(a) * f(b + 1) + f(a - 1) * f(b)) % modulo

    # Thousands of bits, deeper than the default recursion limit
    assert fibonacci(10**1000, 10**9 + 7) == 552179166


if __name__ == '__main__':
    main()
