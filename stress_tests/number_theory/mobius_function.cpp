#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/mobius_function.cpp"
#undef main

int main(){
    /// mu from an Eratosthenes smallest prime factor table, independent of both generators
    int big = 10000000;
    vector<int> spf(big + 1, 0), expected_primes;
    vector<signed char> expected(big + 1, 0);
    for (int i = 2; i <= big; i++){
        if (spf[i]) continue;
        expected_primes.push_back(i);
        for (long long j = i; j <= big; j += i) if (!spf[j]) spf[j] = i;
    }

    expected[1] = 1;
    for (int i = 2; i <= big; i++){
        int p = spf[i], rest = i / p;
        expected[i] = rest % p == 0 ? 0 : -expected[rest];
    }

    MobiusSieve sieve(big);
    assert(sieve.mu == expected && sieve.primes == expected_primes);

    int mid = 1000000;
    assert(mobius_by_divisors(mid) == vector<signed char>(expected.begin(), expected.begin() + mid + 1));

    /// Small sizes, where the sentinel and the bounds are easiest to get wrong
    for (long long it = 0; it < stress::scaled(300); it++){
        int n = it < 3 ? it : stress::rand_int(0, 3000);
        auto prefix = vector<signed char>(expected.begin(), expected.begin() + n + 1);
        assert(MobiusSieve(n).mu == prefix && mobius_by_divisors(n) == prefix);
    }
    return 0;
}
