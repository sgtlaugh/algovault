#include "../common.h"

#define main library_main
#include "../../code_library/misc/divisors.cpp"
#undef main

/// Divisors by testing every candidate up to sqrt(n)
vector<long long> brute(long long n){
    vector<long long> small, large;
    for (long long d = 1; d * d <= n; d++){
        if (n % d) continue;
        small.push_back(d);
        if (d * d != n) large.push_back(n / d);
    }
    small.insert(small.end(), large.rbegin(), large.rend());
    return small;
}

int main(){
    for (long long n = 1; n <= 20000; n++) assert(divisors(n) == brute(n));
    for (long long it = 0; it < stress::scaled(300); it++){
        long long n = it % 3 == 0 ? stress::rand_int(1, 100000000000000LL) : stress::rand_int(1, 1000000000000LL);
        assert(divisors(n) == brute(n));
    }
    /// Highly composite numbers and a prime square at the 1e14 bound
    for (long long n : {963761198400LL, 97821761637600LL, 9999991LL * 9999991LL, 100000000000000LL}) assert(divisors(n) == brute(n));

    /// divisors_from_factors on shuffled factor lists, brute force up to 1e10, structural checks up to 2^63
    const vector<long long> primes = {2, 3, 5, 7, 11, 13, 101, 9973, 1000003, 999999937};
    for (long long it = 0; it < stress::scaled(3000); it++){
        long long limit = it % 2 ? 10000000000LL : LLONG_MAX;
        vector<long long> factors;
        long long n = 1;
        for (int tries = 0; tries < 40; tries++){
            long long p = primes[stress::rand_int(0, it % 3 ? 5 : 9)];
            if (n > limit / p) continue;
            n *= p, factors.push_back(p);
        }
        shuffle(factors.begin(), factors.end(), stress::rng());
        auto got = divisors_from_factors(factors);
        if (limit < LLONG_MAX){
            assert(got == brute(n));
            continue;
        }
        map<long long, int> exponent;
        for (long long p : factors) exponent[p]++;
        long long count = 1;
        for (auto [p, e] : exponent) count *= e + 1;
        assert((long long)got.size() == count && got.front() == 1 && got.back() == n);
        for (size_t i = 0; i < got.size(); i++) assert(n % got[i] == 0 && (i == 0 || got[i - 1] < got[i]));
    }
    assert(divisors_from_factors({}) == vector<long long>{1});
    return 0;
}
