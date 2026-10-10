"""
Berlekamp-Massey

Shortest linear recurrence of a sequence mod an odd prime, then its nth term by Kitamasa

Complexity:
  - berlekamp_massey: O(n^2 + n log mod), n = len(sequence), one modular inverse per update
  - solve_linear_recurrence: O(n^2 + k^2 log(nth_term)), k <= n the recurrence length
"""

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
    if not k:
        return 0

    # Kitamasa: a_N = sum c_i * a_i where sum c_i * x^i = x^N mod the characteristic polynomial
    coefficients = [(mod - r) % mod for r in recurrence]
    poly = [1] + [0] * (k - 1)
    for bit in bin(nth_term)[2:]:
        poly = _multiply_mod(poly, poly, coefficients, mod)
        if bit == '1':
            poly = _multiply_mod(poly, [0, 1], coefficients, mod)

    return sum(c * a for c, a in zip(poly, base_sequence)) % mod


def _multiply_mod(first, second, coefficients, mod):
    """Product of two polynomials reduced by x^k = sum coefficients[j] * x^(k - 1 - j)"""
    k = len(coefficients)
    res = [0] * (len(first) + len(second) - 1)
    for i, f in enumerate(first):
        if f:
            for j, s in enumerate(second):
                res[i + j] += f * s

    for i in range(len(res) - 1, k - 1, -1):
        top = res[i] % mod
        if top:
            for j, c in enumerate(coefficients):
                res[i - 1 - j] += top * c

    return [x % mod for x in res[:k]]


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
