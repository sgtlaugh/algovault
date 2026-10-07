#include "common.h"

#define main library_main
#include "../code_library/binomial_coefficients.cpp"
#undef main

uint64_t random_modulus(){
    if (stress::rand_int(0, 1)) return stress::rand_int(1, 3000);

    /// Up to three coprime prime powers below 1e5, the product reaches ~1e15 where the old CRT overflowed
    uint64_t res = 1;
    vector<int> used;
    for (int parts = stress::rand_int(1, 3); parts; parts--){
        int p;
        do {
            p = stress::rand_int(2, stress::rand_int(0, 1) ? 50 : 100000);
            for (int d = 2; d * d <= p; d++) if (p % d == 0) p = 0;
        } while (!p || count(used.begin(), used.end(), p));
        used.push_back(p);

        uint64_t pq = p;
        while (pq * p < 100000 && stress::rand_int(0, 1)) pq *= p;
        res *= pq;
    }
    return res;
}

/// Rows of Pascal's triangle mod m up to n
void check_pascal(uint64_t m, int n){
    Binomial b(m);
    vector<unsigned __int128> row;
    for (int i = 0; i <= n; i++){
        row.push_back(1 % m);
        for (int k = i - 1; k > 0; k--) row[k] = (row[k] + row[k - 1]) % m;
        for (int k = 0; k <= i + 1; k++) assert((uint64_t)b.binomial(i, k) == (k <= i ? (uint64_t)row[k] : 0));
    }
}

/// Lucas: C(n, k) mod p is the product of C(n_i, k_i) over base p digits
uint64_t lucas(int64_t n, int64_t k, int p){
    vector<vector<uint64_t>> c(p, vector<uint64_t>(p, 0));
    for (int i = 0; i < p; i++){
        c[i][0] = 1;
        for (int j = 1; j <= i; j++) c[i][j] = (c[i - 1][j - 1] + c[i - 1][j]) % p;
    }

    uint64_t res = 1 % p;
    for (; n || k; n /= p, k /= p) res = res * c[n % p][k % p] % p;
    return res;
}

int main(){
    check_pascal(18410739107493357137ULL, 300);  /// 65521 * 65519 * 65497 * 65479, above 2^63

    for (long long it = 0; it < stress::scaled(150); it++){
        uint64_t m = random_modulus();
        check_pascal(m, stress::rand_int(0, 60));

        Binomial b(m);
        for (int q = 0; q < 50; q++){
            int64_t n = stress::rand_int(1, q % 2 ? 1000000 : 1000000000000000000LL), k = stress::rand_int(0, n);
            if (q % 5 == 0) k = stress::rand_int(0, min<int64_t>(n, 20));
            uint64_t a = b.binomial(n, k), l = b.binomial(n - 1, k - 1), r = b.binomial(n - 1, k);
            assert(a == (l >= m - r ? l - (m - r) : l + r));  /// Pascal's rule, exercising fact() across many blocks
            assert(a == (uint64_t)b.binomial(n, n - k));
            assert(a < m || m == 1);
            assert(b.binomial(n, n + 1 + stress::rand_int(0, 5)) == 0);
        }
    }

    const int primes[] = {2, 3, 5, 7, 11, 13, 31, 97};
    for (long long it = 0; it < stress::scaled(400); it++){
        int p = primes[stress::rand_int(0, 7)];
        int64_t n = stress::rand_int(0, 1000000000000000000LL), k = stress::rand_int(0, n);
        assert((uint64_t)Binomial(p).binomial(n, k) == lucas(n, k, p));
    }
    return 0;
}
