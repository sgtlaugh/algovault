/***
 *
 * Modular Square Root
 * Square root modulo a prime with Tonelli-Shanks
 *
 * Complexity: O(log^2 p) expected
 *
 * sqrt_mod(a, p): the smaller of the two x with x^2 = a mod p (the other is p - x), -1 if there is none
 *
 * Notes:
 *   - p must be prime (not checked) and below 2^62, a is any long long and is reduced mod p
 *   - Requires __int128 (64-bit GCC or Clang)
 *   - Multiplicative order, primitive root and k-th roots: modular_roots.cpp, which carries this file
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// BEGIN SHARED modular_sqrt
long long mul_mod(long long a, long long b, long long m){
    return (unsigned __int128)a * b % m;
}

long long pow_mod(long long x, long long e, long long m){
    long long res = 1 % m;
    for (x %= m; e; e >>= 1, x = mul_mod(x, x, m)){
        if (e & 1) res = mul_mod(res, x, m);
    }
    return res;
}

/// Smallest n that is not a q-th power mod p, for a prime q dividing p - 1
long long non_residue(long long q, long long p){
    long long n = 2;
    while (pow_mod(n, (p - 1) / q, p) == 1) n++;
    return n;
}

long long sqrt_mod(long long a, long long p){
    a = (a % p + p) % p;
    if (a == 0 || p == 2) return a;
    if (pow_mod(a, (p - 1) / 2, p) != 1) return -1;

    long long s = p - 1;
    int r = 0;
    for (; s % 2 == 0; s /= 2) r++;

    /// Invariant: x^2 = a * b, and the order of b is 2^m for some m < r
    long long x = pow_mod(a, (s + 1) / 2, p), b = pow_mod(a, s, p), g = pow_mod(non_residue(2, p), s, p);
    while (b != 1){
        int m = 0;
        for (long long t = b; t != 1; t = mul_mod(t, t, p)) m++;

        long long gs = pow_mod(g, 1LL << (r - m - 1), p);
        g = mul_mod(gs, gs, p), x = mul_mod(x, gs, p), b = mul_mod(b, g, p), r = m;
    }

    return min(x, p - x);
}
/// END SHARED modular_sqrt

int main(){
    assert(sqrt_mod(2, 7) == 3);               /// 3^2 = 4^2 = 2 mod 7, the smaller is returned
    assert(sqrt_mod(3, 7) == -1);
    assert(sqrt_mod(0, 11) == 0);
    assert(sqrt_mod(1, 2) == 1);
    assert(sqrt_mod(-1, 13) == 5);             /// 5^2 = 25 = -1 mod 13
    assert(sqrt_mod(-1, 17) == 4);             /// 16 = -1 mod 17, p - 1 = 2^4 exercises the Tonelli-Shanks loop
    assert(sqrt_mod(1000000000000000000LL, (1LL << 62) - 57) == 1000000000);
    return 0;
}
