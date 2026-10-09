#include "../common.h"

#define main library_main
#include "../../code_library/algebra/linear_recurrence.cpp"
#undef main

/// Terms 0 .. count - 1 of f(x) = sum c[i] * f(x - k + i), the coefficient order the constructor takes
vector<int> brute_terms(const vector<int>& c, const vector<int>& start, int count, long long mod){
    int k = c.size();
    vector<int> f = start;
    while ((int)f.size() < count){
        long long v = 0;
        for (int i = 0; i < k; i++) v = (v + (long long)c[i] * f[f.size() - k + i]) % mod;
        f.push_back(v);
    }
    return f;
}

int main(){
    const long long mods[] = {1000000007, 998244353, 1000000, 999999937, 1 << 29, 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23};
    for (long long it = 0; it < stress::scaled(40); it++){
        long long mod = mods[stress::rand_int(0, 5)];
        int k = stress::rand_int(1, it % 4 ? 8 : 40);
        vector<int> c(k), start(k);
        for (auto& x : c) x = stress::rand_int(0, mod - 1);
        for (auto& x : start) x = stress::rand_int(0, mod - 1);
        if (stress::rand_int(0, 2) == 0){  /// zero f(x - k), f(x - k + 1), ... coefficients lower the true order below k
            int zeros = stress::rand_int(1, k);
            fill(c.begin(), c.begin() + zeros, 0);
        }
        auto terms = brute_terms(c, start, 3000, mod);

        /// Recurrence given
        LinearRecurrence given(vector<int>(terms.begin(), terms.begin() + 2 * k), mod, c);
        for (int q = 0; q < 6; q++){
            long long n = stress::rand_int(0, 2999);
            assert(given.nth_term(n) == terms[n]);
        }

        /// Huge n: the k + 1 terms starting at n must satisfy the recurrence themselves
        long long n = stress::rand_int(0, 1000000000000000000LL);
        auto window = given.nth_terms(n, k);
        auto next = given.nth_terms(n + 1, k);
        long long predicted = 0;
        for (int i = 0; i < k; i++) predicted = (predicted + (long long)c[i] * window[i]) % mod;
        assert(next[k - 1] == predicted);
        for (int i = 0; i + 1 < k; i++) assert(next[i] == window[i + 1]);

        /// Recurrence derived from 2k terms (Berlekamp-Massey for primes, Reeds-Sloane otherwise)
        LinearRecurrence derived(vector<int>(terms.begin(), terms.begin() + 2 * k), mod);
        for (int q = 0; q < 6; q++){
            long long m = stress::rand_int(0, 2999);
            assert(derived.nth_term(m) == terms[m]);
        }

        /// Interleaved with a different modulus
        LinearRecurrence other({0, 1, 1, 2}, 7);
        assert(other.nth_term(10) == 55 % 7 && given.nth_term(2 * k + 5) == terms[2 * k + 5] && other.nth_term(11) == 89 % 7);
    }
    return 0;
}
