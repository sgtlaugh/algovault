#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/floor_sum.cpp"
#undef main

__int128 brute_floor(__int128 x, __int128 m){
    __int128 q = x / m;
    return q - (x % m < 0);
}

long long brute_mod(__int128 x, long long m){
    return x - brute_floor(x, m) * m;
}

__int128 brute_floor_sum(long long n, long long m, long long a, long long b){
    __int128 res = 0;
    for (long long i = 0; i < n; i++) res += brute_floor((__int128)a * i + b, m);
    return res;
}

long long brute_min(long long n, long long m, long long a, long long b){
    long long res = m;
    for (long long i = 0; i < n; i++) res = min(res, brute_mod((__int128)a * i + b, m));
    return res;
}

long long brute_count(long long n, long long m, long long a, long long l, long long r){
    long long res = 0;
    for (long long x = 0; x < n; x++){
        long long v = brute_mod((__int128)a * x, m);
        res += l <= v && v <= r;
    }
    return res;
}

/// The residues of a * x repeat with period m, so scanning x < m finds the first hit if any exists
long long brute_first(long long m, long long a, long long l, long long r){
    for (long long x = 0; x < m; x++){
        long long v = brute_mod((__int128)a * x, m);
        if (l <= v && v <= r) return x;
    }
    return -1;
}

/// a * x = v (mod m) is solvable iff g = gcd(a, m) divides v, and the smallest solution is (v / g) * inverse(a / g) mod (m / g)
long long inverse_first(long long m, long long a, long long v){
    __int128 old_r = brute_mod(a, m), r = m, old_s = 1, s = 0;
    while (r){
        __int128 q = old_r / r;
        old_r -= q * r, old_s -= q * s;
        swap(old_r, r), swap(old_s, s);
    }

    if (v % old_r) return -1;
    long long period = m / old_r;
    return brute_mod(v / old_r * old_s, period);
}

long long rand_signed(long long bound){
    return stress::rand_int(-bound, bound);
}

/// Values near 0, near +-m and near the long long limits, where normalization and overflow bugs live
long long rand_edge(long long m){
    __int128 v;
    switch (stress::rand_int(0, 5)){
        case 0: v = stress::rand_int(-3, 3); break;
        case 1: v = (__int128)m + stress::rand_int(-3, 3); break;
        case 2: v = -(__int128)m + stress::rand_int(-3, 3); break;
        case 3: v = LLONG_MAX - stress::rand_int(0, 3); break;
        case 4: v = LLONG_MIN + stress::rand_int(0, 3); break;
        default: v = rand_signed(LLONG_MAX);
    }
    return max<__int128>(LLONG_MIN, min<__int128>(LLONG_MAX, v));
}

void check_all(long long n, long long m, long long a, long long b, long long l, long long r){
    assert(floor_sum(n, m, a, b) == brute_floor_sum(n, m, a, b));
    if (n) assert(min_mod_linear(n, m, a, b) == brute_min(n, m, a, b));
    assert(count_mod_in_range(n, m, a, l, r) == brute_count(n, m, a, max(l, 0LL), min(r, m - 1)));
}

/// A first hit must lie in range with no hit before it; no hit must mean no x in a full period lies in range
void check_first(long long m, long long a, long long l, long long r){
    long long x = first_mod_in_range(m, a, l, r);
    if (x == -1){
        assert(count_mod_in_range(m, m, a, l, r) == 0);
        return;
    }

    long long v = brute_mod((__int128)a * x, m);
    assert(0 <= x && x < m && max(l, 0LL) <= v && v <= min(r, m - 1));
    assert(count_mod_in_range(x, m, a, l, r) == 0);
}

int main(){
    for (long long m = 1; m <= 7; m++){
        for (long long a = -2 * m; a <= 2 * m; a++){
            for (long long b = -2 * m; b <= 2 * m; b++){
                for (long long n = 0; n <= 10; n++) check_all(n, m, a, b, b, b + m / 2);
            }

            for (long long l = -1; l <= m; l++){
                for (long long r = l - 1; r <= m; r++) assert(first_mod_in_range(m, a, l, r) == brute_first(m, a, max(l, 0LL), min(r, m - 1)));
            }
        }
    }

    for (long long it = 0; it < stress::scaled(8000); it++){
        long long m = stress::rand_int(1, it % 2 ? 60 : 3000), n = stress::rand_int(0, it % 3 ? 100 : 4000);
        long long a = rand_signed(3 * m), b = rand_signed(3 * m), l = stress::rand_int(-2, m), r = stress::rand_int(l - 2, m + 2);
        check_all(n, m, a, b, l, r);
        assert(first_mod_in_range(m, a, l, r) == brute_first(m, a, max(l, 0LL), min(r, m - 1)));
    }

    /// Long progressions over moderate moduli drive min_mod_linear through many rounds
    for (long long it = 0; it < stress::scaled(60); it++){
        long long m = stress::rand_int(1, 200000), n = stress::rand_int(1, 200000);
        long long a = stress::rand_int(0, m - 1), b = stress::rand_int(0, m - 1);
        assert(min_mod_linear(n, m, a, b) == brute_min(n, m, a, b));
        assert(floor_sum(n, m, a, b) == brute_floor_sum(n, m, a, b));
    }

    /// Fibonacci ratios maximize the number of Euclid rounds
    long long fib_prev = 1, fib = 1;
    while (fib < 100000){
        long long next = fib + fib_prev;
        fib_prev = fib, fib = next;
        for (long long n : {1LL, fib_prev, fib - 1, fib, 2 * fib}){
            for (long long b : {0LL, 1LL, fib / 2, fib - 1}){
                assert(min_mod_linear(n, fib, fib_prev, b) == brute_min(n, fib, fib_prev, b));
                assert(floor_sum(n, fib, fib_prev, b) == brute_floor_sum(n, fib, fib_prev, b));
                check_first(fib, fib_prev, b, b + stress::rand_int(0, 3));
            }
        }
    }

    /// Full 64-bit moduli and coefficients, short progressions so the loops stay exact references
    for (long long it = 0; it < stress::scaled(4000); it++){
        long long m = it % 4 ? stress::rand_int(1, LLONG_MAX) : LLONG_MAX - stress::rand_int(0, 3);
        long long n = stress::rand_int(0, 300), a = rand_edge(m), b = rand_edge(m);
        long long l = stress::rand_int(-2, m - 1), r = it % 3 ? l + stress::rand_int(0, (m - max(l, 0LL)) / stress::rand_int(1, 1000)) : LLONG_MAX;
        check_all(n, m, a, b, l, r);
        check_first(m, a, l, r);
    }

    /// Known answers at full 64-bit moduli: single residues against the modular inverse, with a shared factor half the time so
    /// unsolvable targets occur, and a = -1 where residue m - x makes the first hit of [l, r] equal m - r
    for (long long it = 0; it < stress::scaled(4000); it++){
        long long g = it % 2 ? stress::rand_int(2, 1000000) : 1;
        long long m = g * stress::rand_int(1, LLONG_MAX / g), a = g * stress::rand_int(0, m / g - 1);
        long long v = it % 4 < 2 ? stress::rand_int(0, m - 1) : g * stress::rand_int(0, m / g - 1);
        assert(first_mod_in_range(m, a, v, v) == inverse_first(m, a, v));
        assert(first_mod_in_range(m, a - m, v, v) == inverse_first(m, a, v));

        if (m == 1) continue;
        long long l = stress::rand_int(1, m - 1), r = stress::rand_int(l, m - 1);
        assert(first_mod_in_range(m, m - 1, l, r) == m - r);
        assert(first_mod_in_range(m, -1, l, r) == m - r);
    }

    /// The header bounds: any a, b when n <= 2^32 (m = 1 makes the sum a closed form), any n when 0 <= a < m
    for (long long it = 0; it < stress::scaled(2000); it++){
        long long n = (1LL << 32) - stress::rand_int(0, 1000), a = rand_edge(1), b = rand_edge(1);
        assert(floor_sum(n, 1, a, b) == (__int128)a * ((__int128)n * (n - 1) / 2) + (__int128)b * n);

        long long big_n = LLONG_MAX - stress::rand_int(0, 1000), k = stress::rand_int(1, 1000000), rb = brute_mod(b, k);
        auto prefix = [&](__int128 x){ __int128 q = x / k; return k * q * (q - 1) / 2 + x % k * q; };
        assert(floor_sum(big_n, k, 1, 0) == prefix(big_n));
        assert(floor_sum(big_n, k, 1, b) == brute_floor(b, k) * big_n + prefix((__int128)big_n + rb) - prefix(rb));
        assert(count_mod_in_range(big_n, k, 1, 0, 0) == ((__int128)big_n + k - 1) / k);

        /// Once n covers a full period the minimum is b mod gcd(a, m)
        long long step = stress::rand_int(0, k - 1);
        assert(min_mod_linear(big_n, k, step, b) == brute_mod(b, __gcd(step, k)));
    }

    return 0;
}
