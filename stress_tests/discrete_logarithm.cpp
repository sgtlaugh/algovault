#include "common.h"

#define main library_main
#include "../code_library/discrete_logarithm.cpp"
#undef main

/// The library's convention: with mod = 1, x = 0 matches only a literal b == 1 and every other b needs x = 1
int brute(int a, int b, int mod){
    if (b == 1) return 0;
    b %= mod;
    long long e = 1;
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
    return 0;
}
