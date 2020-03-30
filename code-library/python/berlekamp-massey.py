from copy import deepcopy


class Matrix:
    def __init__(self, dimension, modulo, diagonal=0):
        """

        :param dimension: size of the matrix, dimension x dimension
        :param modulo: all calculations occur modulo this number
        :param diagonal: diagonal values, matrix[i][i]
        """
        self.n = dimension
        self.modulo = modulo

        self.matrix = [[0 for _ in range(dimension)] for _ in range(dimension)]
        for i in range(dimension):
            self.matrix[i][i] = diagonal

    def multiply(self, multiplicand):
        result = Matrix(self.n, self.modulo)

        for i in range(self.n):
            for j in range(self.n):
                for k in range(self.n):
                    result.matrix[i][j] += self.matrix[i][k] * multiplicand.matrix[k][j]
                result.matrix[i][j] %= self.modulo

        return result

    def exponentiate(self, exponent):
        x = deepcopy(self)
        result = Matrix(self.n, self.modulo, 1)

        while exponent:
            if exponent & 1:
                result = result.multiply(x)

            exponent >>= 1
            x = x.multiply(x)

        return result


def convolution(first, second, modulo):
    assert len(first) == len(second)

    result = 0
    for f, s in zip(first, second):
        result = (result + f * s) % modulo

    return result


def berlekamp_massey(base_sequence, modulo):
    n = len(base_sequence)
    assert n and n % 2 == 0

    U = [int(i == 0) for i in range(n + 1)]
    V = [int(i == 0) for i in range(n + 1)]
    base_sequence = base_sequence[::-1]

    l, m, b, deg = 0, 1, 1, 0
    for i in range(n):
        d = base_sequence[n - i - 1]

        if l:
            d = (d + convolution(V[1:1 + l], base_sequence[n - i:n - i + l], modulo)) % modulo

        if not d:
            m += 1
            continue

        if (l * 2) <= i:
            W = V[:l + 1]

        x = (pow(b, modulo - 2, modulo) * (modulo - d) % modulo + modulo) % modulo
        for j in range(deg + 1):
            V[m + j] = (V[m + j] + x * U[j]) % modulo

        if (l * 2) <= i:
            U, W = W, U
            deg = len(U) - 1
            b, m, l = d, 1, i - l + 1
        else:
            m += 1

    V = V[:l + 1] + [0 for _ in range(max(0, l - len(V) + 1))]
    return V[1:]


def solve_linear_recurrence(base_sequence, nth_term, modulo):
    """

    :param base_sequence: If the recurrence degree is K, len(base_sequence) should be at least 2*K and even
    :param nth_term: nth_term of the recurrence to evaluate
    :param modulo: all calculations will occur modulo this number, needs to be an odd prime
    :return: remainder when the nth_term of the recurrence is divided by modulo
    """

    assert len(base_sequence) % 2 == 0
    base_sequence = [val % modulo for val in base_sequence]

    n = len(base_sequence)
    if nth_term < n:
        return base_sequence[nth_term]

    recurrence = berlekamp_massey(base_sequence, modulo)

    l = len(recurrence)
    ar = Matrix(l, modulo)

    for i in range(l):
        ar.matrix[0][i] = modulo - recurrence[i]
        if i:
            ar.matrix[i][i - 1] = 1

    result = 0
    ar = ar.exponentiate(nth_term - n + 1)
    for i in range(l):
        result += ar.matrix[0][i] * base_sequence[n - i - 1]

    return result % modulo


def main():
    base_sequence = [0, 1, 1, 2, 3, 5, 8, 13]
    print(solve_linear_recurrence(base_sequence, 10, 10 ** 9 + 7))
    print(solve_linear_recurrence(base_sequence, 10 ** 18, 10 ** 9 + 7))


if __name__ == '__main__':
    main()
