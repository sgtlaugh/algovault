#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/du_sieve.cpp"
#undef main

/// Prefix sums of phi and mu from an Eratosthenes smallest prime factor table, independent of the linear sieve
struct Reference{
    vector<long long> phi_sum, mertens;

    Reference(int n) : phi_sum(n + 1, 0), mertens(n + 1, 0){
        vector<int> spf(n + 1, 0);
        for (int i = 2; i <= n; i++){
            if (spf[i]) continue;
            for (int j = i; j <= n; j += i) if (!spf[j]) spf[j] = i;
        }

        vector<long long> phi(n + 1, 0);
        vector<int> mu(n + 1, 0);
        phi[1] = mu[1] = 1;
        for (int i = 2; i <= n; i++){
            int p = spf[i], rest = i / p;
            mu[i] = rest % p == 0 ? 0 : -mu[rest];
            phi[i] = rest % p == 0 ? phi[rest] * p : phi[rest] * (p - 1);
        }

        for (int i = 1; i <= n; i++){
            phi_sum[i] = phi_sum[i - 1] + phi[i];
            mertens[i] = mertens[i - 1] + mu[i];
        }
    }
};

void check(DuSieve& du, const Reference& ref, long long x){
    assert(du.phi_sum(x) == ref.phi_sum[x]);
    assert(du.mertens(x) == ref.mertens[x]);
}

int main(){
    const int big = 10000000;
    Reference ref(big);

    /// Every query on small n, with the smallest tables so nearly everything goes through the recursion
    for (int n = 0; n <= 300; n++){
        for (int table : {1, 2, 3, 7, n}){
            DuSieve du(n, table);
            for (int x = n; x >= 0; x--) check(du, ref, x);
            check(du, ref, 2 * n + 5);
        }

        DuSieve du(n);
        for (int x = 0; x <= n; x++) check(du, ref, x);
    }

    /// The default table at n = 1e7, around its boundary and at random points, in random order on one instance
    DuSieve full(big);
    for (long long x : {full.limit - 1LL, (long long)full.limit, full.limit + 1LL, big - 1LL, (long long)big}) check(full, ref, x);
    for (long long x : {3161LL * 3161 - 1, 3162LL * 3162, 9999991LL, 1LL << 23}) check(full, ref, x);
    for (long long it = 0; it < stress::scaled(300); it++) check(full, ref, stress::rand_int(0, big));

    /// Random n with random small tables, queried at n, at values n / k the recursion visits, and anywhere below n
    for (long long it = 0; it < stress::scaled(20); it++){
        long long n = stress::rand_int(1, big);
        DuSieve du(n, (int)stress::rand_int(1, 2000));
        check(du, ref, n);
        for (int k = 0; k < 20; k++) check(du, ref, n / stress::rand_int(1, n));
        for (int k = 0; k < 4; k++) check(du, ref, stress::rand_int(0, n));
    }

    return 0;
}
