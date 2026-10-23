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

    /// Fewer than 2 terms never confirm a degree, not even the zero polynomial's
    for (int mod : primes){
        assert(find_degree({}, mod) == -1);
        assert(find_degree({0}, mod) == -1);
        assert(find_degree({1}, mod) == -1);
        assert(find_degree({0, 0}, mod) == 0);
        assert(find_degree({-mod, mod}, mod) == 0);
    }

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
        if (stress::rand_int(0, 1)){
            for (auto& t : terms) t -= mod * stress::rand_int(0, 1);
        }

        Lagrange f(terms, mod);
        for (int q = 0; q < 20; q++){
            long long x = q < 5 ? stress::rand_int(-2 * n, 2 * n) : stress::rand_int(-1000000000000000000LL, 1000000000000000000LL);
            assert(f(x) == horner(c, x, mod));
        }

        int expected = zero ? (n >= 2 ? 0 : -1) : (n >= d + 2 ? d : -1);
        assert(find_degree(terms, mod) == expected);
    }

    /// Roots at the first samples and at -1, -2, -3 make wrong lower-degree fits agree with f at small and top residues
    for (long long it = 0; it < stress::scaled(3000); it++){
        int mod = primes[stress::rand_int(0, 10)];
        int d = stress::rand_int(0, min(mod - 2, 12));

        vector<int> c = {(int)stress::rand_int(1, mod - 1)};
        for (int i = 0; i < d; i++){
            int root = stress::rand_int(0, 1) ? stress::rand_int(0, 6) % mod : mod - 1 - stress::rand_int(0, 2) % mod;
            vector<int> next(c.size() + 1, 0);
            for (size_t j = 0; j < c.size(); j++){
                next[j + 1] = (next[j + 1] + c[j]) % mod;
                next[j] = (next[j] + (long long)(mod - root) * c[j]) % mod;
            }
            c = next;
        }

        int n = min<long long>(mod, d + 2 + stress::rand_int(0, 3));
        vector<int> terms(n);
        for (int x = 0; x < n; x++) terms[x] = horner(c, x, mod);
        assert(find_degree(terms, mod) == d);
        assert(Lagrange(terms, mod)(mod - 1) == horner(c, mod - 1, mod));
    }

    return 0;
}
