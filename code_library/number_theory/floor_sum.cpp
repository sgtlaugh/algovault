/***
 *
 * Floor Sum
 * Sums, minima and range counts over the arithmetic progression (a * i + b) mod m, by Euclid-like recursion
 *
 * Complexity: O(log m) per call
 * Requires __int128 (64-bit GCC or Clang)
 *
 * floor_sum(n, m, a, b): sum of floor((a * i + b) / m) for 0 <= i < n, n >= 0, m >= 1, any signed a, b
 *     Exact for any long long a, b when n <= 2^32, and for any n >= 0 when 0 <= a < m
 *     Returns __int128, cast to long long when the sum is known to fit
 * min_mod_linear(n, m, a, b): min of (a * i + b) mod m over 0 <= i < n, n >= 1, any signed a, b
 * count_mod_in_range(n, m, a, l, r): number of 0 <= x < n with l <= (a * x) mod m <= r, any n >= 0
 * first_mod_in_range(m, a, l, r): smallest x >= 0 with l <= (a * x) mod m <= r, or -1 if there is none
 *
 * Arguments are long long; l and r may lie outside [0, m - 1], they are clamped to it
 *
 * Example:
 *   floor_sum(4, 10, 6, 3) = 0 + 0 + 1 + 2 = 3
 *   residues of 6 * x mod 10 are 0, 6, 2, 8, 4, 0, ... so count_mod_in_range(4, 10, 6, 1, 5) = 1 and first_mod_in_range(10, 6, 3, 5) = 4
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

typedef unsigned __int128 u128;

/// Called with 0 <= a, b < m, so every partial sum is at most n * (n - 1) < 2^126 and nothing overflows
u128 floor_sum_unsigned(u128 n, u128 m, u128 a, u128 b){
    u128 res = 0;
    while (true){
        if (a >= m){
            res += n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if (b >= m){
            res += n * (b / m);
            b %= m;
        }

        u128 y_max = a * n + b;
        if (y_max < m) return res;
        n = y_max / m, b = y_max % m;
        swap(m, a);
    }
}

long long floor_div(long long x, long long m){
    return x / m - (x % m < 0);
}

long long mod_floor(long long x, long long m){
    long long r = x % m;
    return r < 0 ? r + m : r;
}

__int128 floor_sum(long long n, long long m, long long a, long long b){
    assert(n >= 0 && m >= 1);

    long long qa = floor_div(a, m), qb = floor_div(b, m);
    long long ra = mod_floor(a, m), rb = mod_floor(b, m);

    __int128 res = (__int128)qa * ((__int128)n * (n - 1) / 2) + (__int128)qb * n;
    return res + (__int128)floor_sum_unsigned(n, m, ra, rb);
}

/// Each round keeps only the candidates that can be minimal: the start and the values right after a wrap when the step is
/// small (2a <= m), the end and the values right before a wrap when it is large; those form a progression mod a or m - a
long long min_mod_linear(long long n, long long m, long long a, long long b){
    assert(n >= 1 && m >= 1);
    a = mod_floor(a, m);
    b = mod_floor(b, m);

    long long res = m;
    while (a){
        if (2 * (__int128)a <= m){
            res = min(res, b);
            long long wraps = ((__int128)a * (n - 1) + b) / m;
            if (!wraps) return res;

            n = wraps;
            b = ((b - m) % a + a) % a;
            long long step = (a - m % a) % a;
            m = a, a = step;
        }
        else{
            long long c = m - a;
            res = min(res, (long long)(((__int128)a * (n - 1) + b) % m));
            if (b >= (__int128)c * n) return res;

            n = ((__int128)c * n - 1 - b) / m + 1;
            a = m % c, b %= c, m = c;
        }
    }

    return min(res, b);
}

long long count_mod_in_range(long long n, long long m, long long a, long long l, long long r){
    assert(n >= 0 && m >= 1);
    l = max(l, 0LL), r = min(r, m - 1);
    if (l > r || !n) return 0;

    a = mod_floor(a, m);
    return floor_sum(n, m, a, m - l) - floor_sum(n, m, a, m - 1 - r);
}

/// Requires 0 <= a < m and 0 <= l <= r < m
long long first_mod_in_range_reduced(long long m, long long a, long long l, long long r){
    if (!a) return l ? -1 : 0;

    __int128 c = ((__int128)l + a - 1) / a;
    if (a * c <= r) return c;

    /// No multiple of a lies in [l, r]; with m = k * a + b, a * x - m * y in [l, r] gives -r <= (b * y) mod a <= -l
    long long b = m % a;
    long long y = first_mod_in_range_reduced(a, b, a - r % a, a - l % a);
    if (y == -1) return -1;
    return ((__int128)b * y + l + a - 1) / a + (__int128)(m / a) * y;
}

long long first_mod_in_range(long long m, long long a, long long l, long long r){
    assert(m >= 1);
    l = max(l, 0LL), r = min(r, m - 1);
    if (l > r) return -1;

    a = mod_floor(a, m);
    return first_mod_in_range_reduced(m, a, l, r);
}

int main(){
    const long long INF = LLONG_MAX;

    assert(floor_sum(4, 10, 6, 3) == 3);
    assert(floor_sum(6, 5, 4, 3) == 13);
    assert(floor_sum(1, 1, 0, 0) == 0);
    assert(floor_sum(31415, 92653, 58979, 32384) == 314095480);
    assert(floor_sum(1000000000, 1000000000, 999999999, 999999999) == 499999999500000000LL);
    assert(floor_sum(0, 7, -5, 3) == 0);
    assert(floor_sum(3, 2, -1, 0) == -2);
    assert(floor_sum(4, 3, -2, 1) == -4);
    assert(floor_sum(5, 4, 3, -9) == -6);
    assert(floor_sum(1LL << 32, 1, LLONG_MIN, LLONG_MIN) == -(((__int128)1 << 126) + ((__int128)1 << 94)));
    assert(floor_sum(INF, 2, 1, 0) == (__int128)((1LL << 62) - 1) * ((1LL << 62) - 1));
    assert(floor_sum(INF, INF, INF - 1, INF - 1) == (__int128)INF * (INF - 1) / 2);
    assert(floor_sum(INF, 1, 0, LLONG_MIN) == (__int128)LLONG_MIN * INF);

    assert(min_mod_linear(1, 10, 6, 3) == 3);
    assert(min_mod_linear(3, 10, 6, 3) == 3);
    assert(min_mod_linear(4, 10, 6, 3) == 1);
    assert(min_mod_linear(2, 7, 3, 5) == 1);
    assert(min_mod_linear(7, 7, 3, 5) == 0);
    assert(min_mod_linear(3, 10, 7, 4) == 1);
    assert(min_mod_linear(8, 10, 7, 4) == 1);
    assert(min_mod_linear(9, 10, 7, 4) == 0);
    assert(min_mod_linear(9, 10, -3, -6) == 0);
    assert(min_mod_linear(5, 10, -3, -6) == 1);
    assert(min_mod_linear(100, 1, 5, 7) == 0);
    assert(min_mod_linear(500000000000000000LL, 1000000000000000000LL, 1, 500000000000000000LL) == 500000000000000000LL);
    assert(min_mod_linear(500000000000000001LL, 1000000000000000000LL, 1, 500000000000000000LL) == 0);
    assert(min_mod_linear(500000000000000000LL, 1000000000000000000LL, -1, 500000000000000000LL) == 1);
    assert(min_mod_linear(INF, 1000000000000000000LL, 600000000000000000LL, 500000000000000007LL) == 100000000000000007LL);

    assert(count_mod_in_range(4, 10, 6, 1, 5) == 1);
    assert(count_mod_in_range(4, 10, 6, 0, 9) == 4);
    assert(count_mod_in_range(4, 10, 6, 6, 8) == 2);
    assert(count_mod_in_range(4, 10, 6, 5, 4) == 0);
    assert(count_mod_in_range(4, 10, 6, -5, 100) == 4);
    assert(count_mod_in_range(0, 10, 6, 0, 9) == 0);
    assert(count_mod_in_range(14, 7, 3, 2, 4) == 6);
    assert(count_mod_in_range(14, 7, -4, 2, 4) == 6);
    assert(count_mod_in_range(1000000000000000000LL, 1000000000000000000LL, 1, 5, 100000000000000000LL) == 99999999999999996LL);
    assert(count_mod_in_range(1000000000000000000LL, 1000000000000000000LL, -1, 0, 0) == 1);
    assert(count_mod_in_range(INF, 2, 1, 1, 1) == INF / 2);

    assert(first_mod_in_range(10, 6, 0, 0) == 0);
    assert(first_mod_in_range(10, 6, 1, 3) == 2);
    assert(first_mod_in_range(10, 6, 3, 5) == 4);
    assert(first_mod_in_range(10, 6, 7, 9) == 3);
    assert(first_mod_in_range(10, 6, 5, 5) == -1);
    assert(first_mod_in_range(10, 6, 9, 3) == -1);
    assert(first_mod_in_range(7, 3, 1, 1) == 5);
    assert(first_mod_in_range(7, -4, 1, 1) == 5);
    assert(first_mod_in_range(1, 5, -3, 3) == 0);
    assert(first_mod_in_range(1000000000000000000LL, 999999999999999999LL, 1, 5) == 999999999999999995LL);
    assert(first_mod_in_range(1000000000000000001LL, 2, 1, 1) == 500000000000000001LL);
    assert(first_mod_in_range(INF, 2, 1, 1) == INF / 2 + 1);

    return 0;
}
