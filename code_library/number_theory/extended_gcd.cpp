/***
 *
 * Extended GCD
 * Extended gcd, modular inverse, and linear diophantine equations (solve and count solutions in a box)
 *
 * Complexity: O(log m) for mod_inverse, O(log min(|a|, |b|)) for the others
 *
 * count_diophantine(a, b, c, x1, x2, y1, y2): number of (x, y) with a * x + b * y = c in [x1, x2] x [y1, y2]
 * diophantine(a, b, c): some (x, y) with a * x + b * y = c, nullopt if there is none
 * extended_gcd(a, b): {g, x, y} with g = gcd(a, b) and a * x + b * y = g, for a, b >= 0
 * mod_inverse(a, m): inverse of a modulo m in [0, m), -1 if gcd(a, m) != 1, for any a and 1 <= m < 2^63
 *
 * Notes:
 *   - count_diophantine: arguments in [-1e9, 1e9], an empty range gives 0, a = b = c = 0 counts the whole box
 *   - diophantine: a, b, c may be negative, |a|, |b|, |c| <= 1e9 keeps x and y within long long
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED extended_gcd
array<long long, 3> extended_gcd(long long a, long long b){
    assert(a >= 0 && b >= 0);

    long long x = 1, y = 0, x1 = 0, y1 = 1;
    while (b){
        long long q = a / b;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(a, b) = make_pair(b, a - q * b);
    }

    return {a, x, y};
}

long long mod_inverse(long long a, long long m){
    assert(m >= 1);

    a %= m;
    if (a < 0) a += m;

    auto [g, x, y] = extended_gcd(a, m);
    if (g != 1) return -1;
    return x < 0 ? x + m : x;
}
/// END SHARED extended_gcd

optional<pair<long long, long long>> diophantine(long long a, long long b, long long c){
    if (a == 0 && b == 0){
        if (c) return nullopt;
        return pair(0LL, 0LL);
    }

    auto [g, x, y] = extended_gcd(llabs(a), llabs(b));
    if (c % g) return nullopt;

    x *= c / g, y *= c / g;
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return pair(x, y);
}

long long count_diophantine(long long a, long long b, long long c, long long x1, long long x2, long long y1, long long y2){
    auto sol = diophantine(a, b, c);
    if (!sol) return 0;
    if (a == 0 && b == 0) return max(0LL, x2 - x1 + 1) * max(0LL, y2 - y1 + 1);

    /// shrinks [lo, hi] to the k with lo_val <= base + k * step <= hi_val
    long long g = gcd(llabs(a), llabs(b)), lo = LLONG_MIN, hi = LLONG_MAX;
    auto restrict_steps = [&](long long base, long long step, long long lo_val, long long hi_val){
        if (step == 0){
            if (base < lo_val || base > hi_val) lo = 1, hi = 0;
            return;
        }

        if (step < 0) step = -step, base = -base, tie(lo_val, hi_val) = make_pair(-hi_val, -lo_val);
        auto floor_div = [&](long long num){ return num / step - (num % step < 0); };
        lo = max(lo, -floor_div(base - lo_val));
        hi = min(hi, floor_div(hi_val - base));
    };
    restrict_steps(sol->first, b / g, x1, x2);
    restrict_steps(sol->second, -a / g, y1, y2);
    return max(0LL, hi - lo + 1);
}

int main(){
    /// Bezout coefficients of 240 and 46, whose gcd is 2
    auto [g, x, y] = extended_gcd(240, 46);
    assert(g == 2);
    assert(240 * x + 46 * y == 2);

    assert(mod_inverse(3, 7) == 5);   /// 3 * 5 = 15 = 1 mod 7
    assert(mod_inverse(-4, 7) == 5);  /// -4 = 3 mod 7
    assert(mod_inverse(4, 6) == -1);  /// gcd(4, 6) = 2

    auto sol = diophantine(3, 5, 7);
    assert(sol);
    assert(3 * sol->first + 5 * sol->second == 7);
    assert(!diophantine(2, 4, 5));  /// gcd(2, 4) does not divide 5

    assert(count_diophantine(2, 3, 12, 0, 10, 0, 10) == 3);  /// (0, 4), (3, 2) and (6, 0)
    assert(count_diophantine(2, -3, 0, -3, 3, -3, 3) == 3);  /// (-3, -2), (0, 0) and (3, 2)
    return 0;
}
