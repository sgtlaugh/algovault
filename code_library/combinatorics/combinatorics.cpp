/***
 *
 * Modular Combinatorics
 * Extended gcd, modular inverse, linear diophantine equations (solve and count in ranges), and factorial based nCr / nPr modulo a prime
 *
 * Complexity:
 *   - O(log m) for extended_gcd, mod_inverse, diophantine and count_diophantine
 *   - O(n) to build Combinatorics(n, mod), O(1) per query afterwards
 *
 * extended_gcd(a, b, x, y): returns g = gcd(a, b) and sets a * x + b * y = g, for a, b >= 0
 * mod_inverse(a, m): inverse of a modulo m for any 1 <= m < 2^63, -1 if gcd(a, m) != 1, a may be negative or >= m
 * diophantine(a, b, c, x, y): solves a * x + b * y = c, returns false if there is no solution
 *     a, b, c may be negative, |a|, |b|, |c| <= 1e9 keeps x and y within long long
 * count_diophantine(a, b, c, x1, x2, y1, y2): number of solutions with x1 <= x <= x2 and y1 <= y <= y2
 *     every argument within [-1e9, 1e9], empty ranges give 0
 *     a = b = 0 and c = 0 gives the full box, (2e9 + 1)^2 at most, still within long long
 *
 * Combinatorics comb(n, mod): tables for 0..n, mod must be a prime with n < mod < 2^31
 * comb.nCr(n, k), comb.nPr(n, k): 0 when k < 0 or k > n
 * comb.inverse(i): inverse of i for 1 <= i <= n, comb.factorial(i), comb.inv_factorial(i)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

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

/// den > 0, restrict_steps flips the sign of step before calling
long long floor_div(long long num, long long den){
    return num / den - (num % den < 0);
}

/// shrinks [lo, hi] to the k with lo_val <= base + k * step <= hi_val
void restrict_steps(long long base, long long step, long long lo_val, long long hi_val, long long& lo, long long& hi){
    if (step == 0){
        if (base < lo_val || base > hi_val) lo = 1, hi = 0;
        return;
    }

    if (step < 0) step = -step, base = -base, tie(lo_val, hi_val) = make_pair(-hi_val, -lo_val);
    lo = max(lo, -floor_div(base - lo_val, step));
    hi = min(hi, floor_div(hi_val - base, step));
}

long long count_diophantine(long long a, long long b, long long c, long long x1, long long x2, long long y1, long long y2){
    long long x, y;
    if (!diophantine(a, b, c, x, y)) return 0;
    if (a == 0 && b == 0) return max(0LL, x2 - x1 + 1) * max(0LL, y2 - y1 + 1);

    long long g = gcd(llabs(a), llabs(b));
    long long lo = LLONG_MIN, hi = LLONG_MAX;
    restrict_steps(x, b / g, x1, x2, lo, hi);
    restrict_steps(y, -a / g, y1, y2, lo, hi);
    return max(0LL, hi - lo + 1);
}

struct Combinatorics{
    long long mod;
    vector<long long> fact, inv_fact, inv;

    Combinatorics(int n, long long mod) : mod(mod), fact(n + 1), inv_fact(n + 1), inv(n + 1){
        assert(n >= 0 && n < mod && mod < (1LL << 31));

        fact[0] = 1;
        for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i % mod;

        if (n >= 1) inv[1] = 1;
        for (int i = 2; i <= n; i++) inv[i] = mod - (mod / i) * inv[mod % i] % mod;  /// needs a prime mod

        inv_fact[0] = 1;
        for (int i = 1; i <= n; i++) inv_fact[i] = inv_fact[i - 1] * inv[i] % mod;
    }

    long long factorial(int i) const{
        return fact[i];
    }

    long long inv_factorial(int i) const{
        return inv_fact[i];
    }

    long long inverse(int i) const{
        assert(i >= 1);
        return inv[i];
    }

    long long nCr(int n, int k) const{
        if (k < 0 || k > n) return 0;
        return fact[n] * inv_fact[k] % mod * inv_fact[n - k] % mod;
    }

    long long nPr(int n, int k) const{
        if (k < 0 || k > n) return 0;
        return fact[n] * inv_fact[n - k] % mod;
    }
};

int main(){
    long long x, y;
    assert(extended_gcd(240, 46, x, y) == 2 && 240 * x + 46 * y == 2);
    assert(extended_gcd(7, 0, x, y) == 7 && 7 * x == 7);
    assert(extended_gcd(0, 0, x, y) == 0);

    assert(mod_inverse(3, 7) == 5);
    assert(mod_inverse(-4, 7) == 5);
    assert(mod_inverse(10, 7) == 5);
    assert(mod_inverse(4, 6) == -1);
    assert(mod_inverse(5, 1) == 0);
    assert(mod_inverse(2, 1000000007) == 500000004);

    assert(diophantine(3, 5, 7, x, y) && 3 * x + 5 * y == 7);
    assert(diophantine(-6, 4, 2, x, y) && -6 * x + 4 * y == 2);
    assert(!diophantine(2, 4, 5, x, y));
    assert(diophantine(0, 0, 0, x, y));
    assert(!diophantine(0, 0, 3, x, y));

    assert(count_diophantine(2, 3, 12, 0, 10, 0, 10) == 3);
    assert(count_diophantine(-2, -3, -12, 0, 10, 0, 10) == 3);
    assert(count_diophantine(1, 1, 5, -5, 10, 2, 4) == 3);
    assert(count_diophantine(-10, -8, -80, -100, 100, -90, 90) == 37);
    assert(count_diophantine(2, 3, 4, 1, 7, 0, 8) == 1);
    assert(count_diophantine(-2, -3, -6, -2, 5, -10, 5) == 2);
    assert(count_diophantine(6, -4, 2, -20, 20, -20, 20) == 14);
    assert(count_diophantine(0, 3, 6, -2, 2, 0, 5) == 5);
    assert(count_diophantine(4, 0, -8, -3, 3, -1, 6) == 8);
    assert(count_diophantine(0, 0, 0, 1, 3, -2, 1) == 12);
    assert(count_diophantine(0, 0, 5, 1, 3, -2, 1) == 0);
    assert(count_diophantine(2, 4, 5, -50, 50, -50, 50) == 0);
    assert(count_diophantine(2, 3, 12, 5, 4, 0, 10) == 0);
    assert(count_diophantine(0, 0, 0, 0, 5, 3, 2) == 0);
    assert(count_diophantine(0, 0, 0, 5, 3, 2, 0) == 0);
    assert(count_diophantine(1, 1, 0, -1000000000, 1000000000, -1000000000, 1000000000) == 2000000001);
    assert(count_diophantine(1000000000, 999999999, 1000000000, -1000000000, 1000000000, -1000000000, 1000000000) == 3);
    assert(count_diophantine(0, 0, 0, -1000000000, 1000000000, -1000000000, 1000000000) == 4000000004000000001LL);

    Combinatorics comb(10, 1000000007);
    assert(comb.nCr(5, 2) == 10);
    assert(comb.nCr(10, 3) == 120);
    assert(comb.nCr(0, 0) == 1);
    assert(comb.nCr(5, 6) == 0);
    assert(comb.nCr(5, -1) == 0);
    assert(comb.nPr(5, 2) == 20);
    assert(comb.nPr(4, 4) == 24);
    assert(comb.factorial(10) == 3628800);
    assert(comb.inverse(2) == 500000004);

    Combinatorics small(6, 7);
    assert(small.nCr(6, 3) == 6);
    assert(small.inverse(3) == 5);
    assert(small.inv_factorial(6) == 6);

    return 0;
}
