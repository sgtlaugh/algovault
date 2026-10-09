"""
Sieve of Eratosthenes

sieve(n): bytearray where is_prime[i] is 1 when i is prime, for 0 <= i <= n
primes_up_to(n): list of the primes <= n

Complexity: O(n log log n), n + 1 bytes for the result plus a temporary of up to n / 2 bytes while marking
Slice assignment marks all multiples of a prime in one C level step, far faster than a Python loop
"""


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


def primes_up_to(n):
    return [i for i, flag in enumerate(sieve(n)) if flag]


def main():
    assert primes_up_to(30) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert primes_up_to(1) == [] and primes_up_to(0) == [] and primes_up_to(2) == [2] and primes_up_to(3) == [2, 3]
    assert len(primes_up_to(1000000)) == 78498
    assert sieve(97)[97] == 1 and sieve(91)[91] == 0 and sieve(4)[4] == 0


if __name__ == '__main__':
    main()
