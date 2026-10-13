"""
Derangements
Ways to seat n people in m seats with person x never in seat x, and the probability a permutation is a derangement

Complexity: O(n) big integer operations on O(n log m) bit values for derangement(n, m), O(1) for
derangement_probability(n) since it returns 1 / e for n > 20
"""

from math import exp, factorial
from decimal import Decimal


def derangement(n, m):
    """
    number of ways to assign n people to m seats subject to the following condition:
    a seat can allocate at most one people and no people numbered x can be assigned to the seat numbered x
    the people are numbered from 1 to n and seats are numbered from 1 to m
    if m=n, then it is the standard derangement problem

    :param n: number of people
    :param m: number of seats
    :return:  number of ways to assign n seats from m seats to n people satisfying the condition above
    """

    if n > m:
        return 0

    x, y = 1, m - n
    for k in range(2, n + 1):
        z = (m - n) * y + (k - 1) * (x + y)
        x, y = y, z

    return x if n == 0 else y


def derangement_probability(n):
    """
    probability that if a random permutation from 1 to n is chosen, it is a derangement
    limit of ratio of derangement to permutation as n approaches infinity is 1 / e
    """
    if n > 20:
        return exp(-1)

    x = derangement(n, n)
    y = factorial(n)
    return float(Decimal(x) / Decimal(y))


def main():
    import sys

    data = sys.stdin.read().split()
    out = []
    for token in data[1:int(data[0]) + 1]:
        n = int(token)
        out.append('%d/%d' % (derangement(n, n), factorial(n)))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
