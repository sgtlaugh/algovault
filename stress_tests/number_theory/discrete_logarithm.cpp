#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/discrete_logarithm.cpp"
#undef main

int residue(long long v, int mod){
    return (v % mod + mod) % mod;
}

int pow_mod(long long x, int n, int mod){
    long long res = 1 % mod;
    for (x = residue(x, mod); n; n >>= 1, x = x * x % mod){
        if (n & 1) res = res * x % mod;
    }
    return res;
}

int brute(int a, int b, int mod){
    a = residue(a, mod), b = residue(b, mod);
    long long e = 1 % mod;
    for (int x = 0; x <= mod; x++, e = e * a % mod){
        if (e == b) return x;
    }
    return -1;
}

int main(){
    for (long long it = 0; it < stress::scaled(20000); it++){
        int mod = stress::rand_int(1, it % 10 ? 100 : 5000);
        int a = stress::rand_int(it % 5 ? 0 : -3 * mod, it % 3 ? mod - 1 : 3 * mod);
        int b = it % 2 ? pow_mod(a, stress::rand_int(0, 2 * mod), mod) : stress::rand_int(-3 * mod, 3 * mod);
        if (it % 4 == 1) b -= mod * stress::rand_int(0, 3);
        assert(discrete_log(a, b, mod) == brute(a, b, mod));
    }

    /// Moduli too large for brute force: b = a^k, so a solution no larger than k must come back
    for (long long it = 0; it < stress::scaled(30); it++){
        int mod = stress::rand_int(it % 2 ? 1000000 : 2000000000, INT_MAX), a = stress::rand_int(it % 3 ? 0 : INT_MIN, mod - 1), k = stress::rand_int(0, mod);
        int b = pow_mod(a, k, mod), x = discrete_log(a, b, mod);
        assert(0 <= x && x <= k && pow_mod(a, x, mod) == b);
    }

    /// INT_MIN for both a and b, whose reduction would overflow if negated first
    assert(discrete_log(INT_MIN, INT_MIN, INT_MAX) == 1);

    /// Baby steps b * 3^i mod 2^30 share their low bits, the worst case for an identity hash with masked buckets
    auto start = chrono::steady_clock::now();
    for (int shift = 11; shift <= 14; shift++){
        for (int odd = 5; odd <= 13; odd += 2) assert(discrete_log(3, odd << shift, 1 << 30) == -1);
    }
    assert(chrono::steady_clock::now() - start < chrono::seconds(10));

    return 0;
}
