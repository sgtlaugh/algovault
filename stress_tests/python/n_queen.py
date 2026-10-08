from itertools import permutations

import stress  # noqa: F401, makes code_library/python importable

from n_queen import count_ways

OEIS = [1, 1, 0, 0, 2, 10, 4, 40, 92, 352, 724, 2680, 14200]  # A000170


def brute(n):
    """One queen per row and column is a permutation, it only remains to rule out shared diagonals"""
    return sum(len({r + c for r, c in enumerate(p)}) == n == len({r - c for r, c in enumerate(p)}) for p in permutations(range(n)))


def main():
    for n in range(8):
        assert brute(n) == OEIS[n]
    for n in range(len(OEIS)):
        assert count_ways(n) == OEIS[n]


if __name__ == '__main__':
    main()
