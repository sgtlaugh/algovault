from copy import deepcopy


class Matrix:
    def __init__(self, n, mod, diagonal=0):
        """

        :param n: size of the matrix, n x n
        :param mod: all calculations occur modulo this number
        :param diagonal: diagonal values, mat[i][i]
        """
        self.n = n
        self.mod = mod

        self.mat = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(n):
            self.mat[i][i] = diagonal

    def multiply(self, multiplicand):
        res = Matrix(self.n, self.mod)

        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    res.mat[i][j] += self.mat[i][k] * multiplicand.mat[k][j]
                res.mat[i][j] %= self.mod

        return res

    def power(self, exponent):
        mat = deepcopy(self)
        res = Matrix(self.n, self.mod, 1)

        while exponent:
            if exponent & 1:
                res = res.multiply(mat)

            exponent >>= 1
            mat = mat.multiply(mat)

        return res


def convolution(first, second, mod):
    res = 0
    for f, s in zip(first, second):
        res = (res + f * s) % mod
    return res


def berlekamp_massey(sequence, mod):
    n = len(sequence)
    assert n and n % 2 == 0

    u_vals = [int(i == 0) for i in range(n + 1)]
    v_vals = [int(i == 0) for i in range(n + 1)]
    sequence = sequence[::-1]

    l, m, b, deg = 0, 1, 1, 0
    for i in range(n):
        d = sequence[n - i - 1]

        if l:
            d = (d + convolution(v_vals[1:1 + l], sequence[n - i:n - i + l], mod)) % mod

        if not d:
            m += 1
            continue

        if (l * 2) <= i:
            w_vals = v_vals[:l + 1]

        x = (pow(b, mod - 2, mod) * (mod - d) % mod + mod) % mod
        for j in range(deg + 1):
            v_vals[m + j] = (v_vals[m + j] + x * u_vals[j]) % mod

        if (l * 2) <= i:
            u_vals, w_vals = w_vals, u_vals
            deg = len(u_vals) - 1
            b, m, l = d, 1, i - l + 1
        else:
            m += 1

    v_vals = v_vals[:l + 1] + [0 for _ in range(max(0, l - len(v_vals) + 1))]
    return v_vals[1:]


def solve_linear_recurrence(base_sequence, nth_term, mod):
    """

    :param base_sequence: If the recurrence degree is k, len(base_sequence) should be at least 2*k and even
    :param nth_term: nth_term of the recurrence to evaluate
    :param mod: all calculations will occur mod this number, needs to be an odd prime
    :return: remainder when the nth_term of the recurrence is divided by mod
    """

    base_sequence = [val % mod for val in base_sequence]

    n = len(base_sequence)
    if n % 2 != 0:
        raise ValueError('Length of base_sequence must be even')

    if nth_term < n:
        return base_sequence[nth_term]

    recurrence = berlekamp_massey(base_sequence, mod)

    k = len(recurrence)
    ar = Matrix(k, mod)

    for i in range(k):
        ar.mat[0][i] = mod - recurrence[i]
        ar.mat[i][i - 1] = int(i > 0)

    result = 0
    ar = ar.power(nth_term - n + 1)
    for i in range(k):
        result += ar.mat[0][i] * base_sequence[n - i - 1]

    return result % mod


def main():
    mod = 10**9 + 7
    base_sequence = [0, 1, 1, 2, 3, 5, 8, 13]  # fibonacci sequence
    
    assert solve_linear_recurrence(base_sequence, 0, mod) == 0
    assert solve_linear_recurrence(base_sequence, 1, mod) == 1
    assert solve_linear_recurrence(base_sequence, 10, mod) == 55
    assert solve_linear_recurrence(base_sequence, 11, mod) == 89
    assert solve_linear_recurrence(base_sequence, 10**18, mod) == 209783453


if __name__ == '__main__':
    main()
