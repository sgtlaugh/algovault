#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/linear_sieve.cpp"
#undef main

/// Trial division factorization, independent of the sieve
vector<pair<long long, int>> trial(long long x){
    vector<pair<long long, int>> res;
    for (long long p = 2; p * p <= x; p++){
        int e = 0;
        for (; x % p == 0; x /= p) e++;
        if (e) res.push_back({p, e});
    }

    if (x > 1) res.push_back({x, 1});
    return res;
}

void check_value(const LinearSieve& sieve, const vector<long long>& divisors, const vector<long long>& sigma, const vector<long long>& identity, int x){
    auto f = trial(x);
    long long phi = x, d = 1, s = 1;
    int mu = 1;
    vector<int> flat;

    for (auto [p, e] : f){
        phi = phi / p * (p - 1);
        d *= e + 1;
        long long pk = 1, sum = 1;
        for (int i = 0; i < e; i++) pk *= p, sum += pk;
        s *= sum;
        mu = e > 1 ? 0 : -mu;
        for (int i = 0; i < e; i++) flat.push_back(p);
    }

    assert(sieve.phi[x] == phi && sieve.mu[x] == mu);
    assert(divisors[x] == d && sigma[x] == s && identity[x] == x);
    assert(sieve.factorize(x) == flat);
    if (x >= 2) assert(sieve.spf[x] == flat[0] && sieve.spf_exponent[x] == f[0].second);
}

/// g(p^k) = p^k gives g[x] = x, and the assert checks that pk really is p^k
vector<long long> identity_of(const LinearSieve& sieve){
    return sieve.multiplicative([](long long p, int k, long long pk){
        long long power = 1;
        for (int i = 0; i < k; i++) power *= p;
        assert(pk == power);
        return pk;
    });
}

int main(){
    int small = 5000;
    LinearSieve a(small);
    auto da = a.multiplicative([](long long, int k, long long){ return (long long)k + 1; });
    auto sa = a.multiplicative([](long long p, int, long long pk){ return (pk * p - 1) / (p - 1); });
    auto ia = identity_of(a);
    for (int x = 1; x <= small; x++) check_value(a, da, sa, ia, x);

    vector<int> brute_primes;
    for (int x = 2; x <= small; x++){
        if (trial(x).size() == 1 && trial(x)[0].second == 1) brute_primes.push_back(x);
    }
    assert(a.primes == brute_primes);

    /// phi via gcd counting on a smaller range, a definition independent of factorization
    for (int x = 1; x <= 600; x++){
        int coprime = 0;
        for (int y = 1; y <= x; y++) coprime += __gcd(x, y) == 1;
        assert(a.phi[x] == coprime);
    }

    int big = 10000000;
    LinearSieve b(big);
    auto db = b.multiplicative([](long long, int k, long long){ return (long long)k + 1; });
    auto sb = b.multiplicative([](long long p, int, long long pk){ return (pk * p - 1) / (p - 1); });
    auto ib = identity_of(b);
    assert(b.primes.size() == 664579);
    for (long long it = 0; it < stress::scaled(3000); it++) check_value(b, db, sb, ib, stress::rand_int(1, big));
    for (int x : {big, big - 1, 9999991, 1 << 23, 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3}) check_value(b, db, sb, ib, x);

    for (int n = 0; n <= 3; n++){
        LinearSieve t(n);
        assert((int)t.phi.size() == n + 1);
    }

    return 0;
}
