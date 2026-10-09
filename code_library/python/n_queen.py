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

    if not n:
        return 1

    # Mirror symmetry: first row columns in the left half count twice, the middle column of an odd board once
    place = lambda x: backtrack(n - 1, x, x << 1, x >> 1)
    ways = 2 * sum(place(1 << i) for i in range(n // 2))
    return ways + place(1 << (n // 2)) if n & 1 else ways


def main():
    assert count_ways(8) == 92
    assert count_ways(13) == 73712


if __name__ == '__main__':
    main()
