"""
Sieve of Eratosthenes

primes_up_to(n): list of the primes <= n
sieve(n): bytearray where is_prime[i] is 1 when i is prime, for 0 <= i <= n

Complexity: O(n log log n), n + 1 bytes for the result plus a temporary of up to n / 2 bytes while marking
Slice assignment marks all multiples of a prime in one C level step, far faster than a Python loop
"""

from itertools import compress


def primes_up_to(n):
    return list(compress(range(n + 1), sieve(n)))


def sieve(n):
    is_prime = bytearray([1]) * (n + 1)
    is_prime[:2] = bytearray(min(2, n + 1))
    if n >= 4:
        is_prime[4::2] = bytearray(len(range(4, n + 1, 2)))

    i = 3
    while i * i <= n:
        if is_prime[i]:
            is_prime[i * i::2 * i] = bytearray(len(range(i * i, n + 1, 2 * i)))
        i += 2

    return is_prime


def main():
    import sys

    primes = primes_up_to(20000000)
    twins = [p for p, q in zip(primes, primes[1:]) if q - p == 2][:100000]

    out = []
    for token in sys.stdin.buffer.read().split():
        p = twins[int(token) - 1]
        out.append('(%d, %d)' % (p, p + 2))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
