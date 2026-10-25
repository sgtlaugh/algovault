/***
 *
 * Extended GCD
 * Extended gcd, modular inverse, and linear diophantine equations (solve and count solutions in a box)
 *
 * Complexity: O(log m) for extended_gcd, mod_inverse, diophantine and count_diophantine
 *
 * count_diophantine(a, b, c, x1, x2, y1, y2): number of solutions with x1 <= x <= x2 and y1 <= y <= y2
 *     every argument within [-1e9, 1e9], empty ranges give 0
 *     a = b = 0 and c = 0 gives the full box, (2e9 + 1)^2 at most, still within long long
 * diophantine(a, b, c, x, y): solves a * x + b * y = c, returns false if there is no solution
 *     a, b, c may be negative, |a|, |b|, |c| <= 1e9 keeps x and y within long long
 * extended_gcd(a, b, x, y): returns g = gcd(a, b) and sets a * x + b * y = g, for a, b >= 0
 * mod_inverse(a, m): inverse of a modulo m for any 1 <= m < 2^63, -1 if gcd(a, m) != 1, a may be negative or >= m
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED extended_gcd
long long extended_gcd(long long a, long long b, long long& x, long long& y){
    assert(a >= 0 && b >= 0);

    long long x1 = 0, y1 = 1;
    x = 1, y = 0;
    while (b){
        long long q = a / b;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(a, b) = make_pair(b, a - q * b);
    }
    return a;
}

long long mod_inverse(long long a, long long m){
    assert(m >= 1);

    a %= m;
    if (a < 0) a += m;

    long long x, y;
    if (extended_gcd(a, m, x, y) != 1) return -1;
    x %= m;
    return x < 0 ? x + m : x;
}
/// END SHARED extended_gcd

bool diophantine(long long a, long long b, long long c, long long& x, long long& y){
    if (a == 0 && b == 0){
        x = y = 0;
        return c == 0;
    }

    long long g = extended_gcd(llabs(a), llabs(b), x, y);
    if (c % g) return false;

    x *= c / g, y *= c / g;
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return true;
}

long long count_diophantine(long long a, long long b, long long c, long long x1, long long x2, long long y1, long long y2){
    long long x, y;

    if (!diophantine(a, b, c, x, y)) return 0;
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
    restrict_steps(x, b / g, x1, x2);
    restrict_steps(y, -a / g, y1, y2);
    return max(0LL, hi - lo + 1);
}

int main(){
    long long x, y;
    assert(extended_gcd(240, 46, x, y) == 2);
    assert(240 * x + 46 * y == 2);
    assert(extended_gcd(7, 0, x, y) == 7 && x == 1 && y == 0);

    assert(mod_inverse(3, 7) == 5);   /// 3 * 5 = 15 = 1 mod 7
    assert(mod_inverse(-4, 7) == 5);  /// -4 = 3 mod 7
    assert(mod_inverse(4, 6) == -1);  /// gcd(4, 6) = 2

    assert(diophantine(3, 5, 7, x, y));
    assert(3 * x + 5 * y == 7);
    assert(!diophantine(2, 4, 5, x, y));                       /// gcd(2, 4) does not divide 5
    assert(count_diophantine(2, 3, 12, 0, 10, 0, 10) == 3);    /// (0, 4), (3, 2) and (6, 0)
    assert(count_diophantine(2, -3, 0, -3, 3, -3, 3) == 3);    /// (-3, -2), (0, 0) and (3, 2)
    assert(count_diophantine(0, 0, 0, 1, 2, 1, 3) == 6);
    return 0;
}
