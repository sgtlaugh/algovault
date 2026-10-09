/***
 *
 * Modular Combinatorics
 * Extended gcd, modular inverse, linear diophantine equations, and factorial based nCr / nPr modulo a prime
 *
 * Complexity:
 *   - O(log m) for extended_gcd, mod_inverse and diophantine
 *   - O(n) to build Combinatorics(n, mod), O(1) per query afterwards
 *
 * extended_gcd(a, b, x, y): returns g = gcd(a, b) and sets a * x + b * y = g, for a, b >= 0
 * mod_inverse(a, m): inverse of a modulo m for any 1 <= m < 2^63, -1 if gcd(a, m) != 1, a may be negative or >= m
 * diophantine(a, b, c, x, y): solves a * x + b * y = c, returns false if there is no solution
 *     a, b, c may be negative, |a|, |b|, |c| <= 1e9 keeps x and y within long long
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
    if (extended_gcd(a, m, x, y) != 1) return m == 1 ? 0 : -1;
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
