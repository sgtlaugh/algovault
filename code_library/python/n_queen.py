def count_ways(n):
    """
    Parameters
        ----------
        n : Dimension of chessboard (n x n)

    Returns
        ----------
        The number of ways to place n queens on the chessboard,
        such that no pair of queens attack each other
    """

    def backtrack(i, c, l, r):
        if not i:
            return 1

        ways, mask = 0, ((1 << n) - 1) & ~(l | r | c)
        while mask:
            x = -mask & mask
            mask ^= x
            ways += backtrack(i - 1, c | x, (l | x) << 1, (r | x) >> 1)

        return ways

    return backtrack(n, 0, 0, 0)


def main():
    assert count_ways(8) == 92
    assert count_ways(13) == 73712


if __name__ == '__main__':
    main()
