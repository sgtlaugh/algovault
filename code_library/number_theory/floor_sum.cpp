/***
 *
 * Floor Sum
 * Sums, minima and range counts over the arithmetic progression (a * i + b) mod m, by Euclid-like reduction
 *
 * Complexity: O(log m) per call
 *
 * count_mod_in_range(n, m, a, l, r): number of 0 <= x < n with l <= (a * x) mod m <= r, any n >= 0
 * first_mod_in_range(m, a, l, r): smallest x >= 0 with l <= (a * x) mod m <= r, or -1 if there is none
 * floor_sum(n, m, a, b): sum of floor((a * i + b) / m) for 0 <= i < n, n >= 0, m >= 1, any signed a, b
 *     Exact for any long long a, b when n <= 2^32, and for any n >= 0 when 0 <= a < m
 *     Returns __int128, cast to long long when the sum is known to fit
 * min_mod_linear(n, m, a, b): min of (a * i + b) mod m over 0 <= i < n, n >= 1, any signed a, b
 *
 * Notes:
 *   - Arguments are long long, l and r may lie outside [0, m - 1], they are clamped to it
 *   - Count of 0 <= x < n, y >= 0 with a x + b y <= c is floor_sum(n, b, -a, c + b) if a >= 0 and a (n - 1) <= c + b
 *   - Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

typedef unsigned __int128 u128;

/// Called with 0 <= a, b < m, so a * n + b < m * (n + 1) and the sum stays below n^2: 64 bits suffice for n, m < 2^32
template<typename U>
U floor_sum_unsigned(U n, U m, U a, U b){
    U res = 0;

    while (true){
        if (a >= m){
            res += n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if (b >= m){
            res += n * (b / m);
            b %= m;
        }

        U y_max = a * n + b;
        if (y_max < m) return res;
        n = y_max / m, b = y_max % m;
        swap(m, a);
    }
}

long long floor_sum_reduce_mod(long long x, long long m){
    x %= m;
    return x < 0 ? x + m : x;
}

__int128 floor_sum(long long n, long long m, long long a, long long b){
    assert(n >= 0 && m >= 1);

    long long ra = floor_sum_reduce_mod(a, m), rb = floor_sum_reduce_mod(b, m);
    long long qa = a / m - (a % m < 0), qb = b / m - (b % m < 0);
    __int128 res = (__int128)qa * ((__int128)n * (n - 1) / 2) + (__int128)qb * n;

    /// The 64-bit path is 46% faster for n, m <= 1e9 and cannot overflow below 2^32, see floor_sum_unsigned
    if (n < (1LL << 32) && m < (1LL << 32)) return res + floor_sum_unsigned<unsigned long long>(n, m, ra, rb);
    return res + (__int128)floor_sum_unsigned<u128>(n, m, ra, rb);
}

long long count_mod_in_range(long long n, long long m, long long a, long long l, long long r){
    assert(n >= 0 && m >= 1);
    l = max(l, 0LL), r = min(r, m - 1);
    if (l > r || !n) return 0;

    a = floor_sum_reduce_mod(a, m);
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

    a = floor_sum_reduce_mod(a, m);
    return first_mod_in_range_reduced(m, a, l, r);
}

/// Each round keeps only the candidates that can be minimal, a progression mod a when 2a <= m, else mod m - a
long long min_mod_linear(long long n, long long m, long long a, long long b){
    assert(n >= 1 && m >= 1);
    a = floor_sum_reduce_mod(a, m);
    b = floor_sum_reduce_mod(b, m);

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

int main(){
    /// (6 * i + 3) / 10 for i = 0..3 is 0, 0, 1, 2
    assert(floor_sum(4, 10, 6, 3) == 3);
    assert(floor_sum(3, 2, -1, 0) == -2);  /// 0 + floor(-1 / 2) + floor(-2 / 2), negative a rounds down

    /// (6 * i + 3) mod 10 for i = 0..3 is 3, 9, 5, 1
    assert(min_mod_linear(3, 10, 6, 3) == 3);
    assert(min_mod_linear(4, 10, 6, 3) == 1);

    /// 6 * x mod 10 for x = 0, 1, 2, ... is 0, 6, 2, 8, 4, 0, ...
    assert(count_mod_in_range(4, 10, 6, 6, 8) == 2);     /// 6 and 8
    assert(count_mod_in_range(4, 10, 6, -5, 100) == 4);  /// l and r are clamped to [0, 9]
    assert(first_mod_in_range(10, 6, 3, 5) == 4);
    assert(first_mod_in_range(10, 6, 5, 5) == -1);       /// residues are always even
    return 0;
}
