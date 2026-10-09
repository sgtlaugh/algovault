#include "../common.h"

#define main library_main
#include "../../code_library/algebra/lagrange_interpolation.cpp"
#undef main

int horner(const vector<int>& c, long long x, int mod){
    x = (x % mod + mod) % mod;
    long long v = 0;
    for (int i = c.size() - 1; i >= 0; i--) v = (v * x + c[i]) % mod;
    return v;
}

int main(){
    const int primes[] = {2, 3, 5, 7, 13, 101, 1009, 65537, 998244353, 1000000007, 2147483647};

    for (long long it = 0; it < stress::scaled(1500); it++){
        int mod = primes[stress::rand_int(0, 10)];
        int d = stress::rand_int(0, min(mod - 1, it % 10 ? 30 : 400));
        int n = min<long long>(mod, d + 1 + stress::rand_int(0, 20));

        vector<int> c(d + 1);
        for (auto& x : c) x = stress::rand_int(0, mod - 1);
        if (!c[d]) c[d] = 1;
        bool zero = stress::rand_int(0, 20) == 0;
        if (zero) c.assign(d + 1, 0);

        vector<int> terms(n);
        for (int x = 0; x < n; x++) terms[x] = horner(c, x, mod);

        auto lagrange = Lagrange(terms, mod);
        for (int q = 0; q < 20; q++){
            long long x = q < 5 ? stress::rand_int(-2 * n, 2 * n) : stress::rand_int(-1000000000000000000LL, 1000000000000000000LL);
            assert(lagrange.interpolate(x) == horner(c, x, mod));
        }

        int l = stress::rand_int(-50, 50);
        auto window = get_terms(terms, mod, l, l + 10);
        for (int i = 0; i <= 10; i++) assert(window[i] == horner(c, l + i, mod));

        int expected = zero ? (n >= 2 ? 0 : -1) : (n >= d + 2 ? d : -1);
        assert(find_degree(terms, mod) == expected);
    }

    return 0;
}
