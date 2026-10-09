#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/discrete_logarithm.cpp"
#undef main

int brute(int a, int b, int mod){
    b %= mod;
    long long e = 1 % mod;
    for (int x = 0; x <= mod; x++, e = e * a % mod){
        if (e == b) return x;
    }
    return -1;
}

int main(){
    for (long long it = 0; it < stress::scaled(20000); it++){
        int mod = stress::rand_int(1, it % 10 ? 100 : 5000);
        int a = stress::rand_int(0, it % 3 ? mod - 1 : 3 * mod);
        int b = it % 2 ? expo(a, stress::rand_int(0, 2 * mod), mod) : stress::rand_int(0, 3 * mod);
        assert(discrete_log(a, b, mod) == brute(a, b, mod));
    }

    /// Moduli too large for brute force: b = a^k, so a solution no larger than k must come back
    for (long long it = 0; it < stress::scaled(30); it++){
        int mod = stress::rand_int(it % 2 ? 1000000 : 2000000000, INT_MAX), a = stress::rand_int(0, mod - 1), k = stress::rand_int(0, mod);
        int b = expo(a, k, mod), x = discrete_log(a, b, mod);
        assert(0 <= x && x <= k && expo(a, x, mod) == b);
    }

    /// Baby steps b * 3^i mod 2^30 share their low bits, an identity hash with masked buckets takes 90 s here instead of 1 s
    auto start = chrono::steady_clock::now();
    for (int shift = 11; shift <= 14; shift++){
        for (int odd = 5; odd <= 13; odd += 2) assert(discrete_log(3, odd << shift, 1 << 30) == -1);
    }
    assert(chrono::steady_clock::now() - start < chrono::seconds(10));
    return 0;
}
