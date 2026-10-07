#include "common.h"

#define main library_main
#include "../code_library/fast_sieve.cpp"
#undef main

/// Primes in [lo, hi] by a plain segmented sieve over the primes below 2^16
vector<uint32_t> segment(uint32_t lo, uint32_t hi, const vector<uint32_t>& small){
    vector<char> composite(hi - lo + 1, 0);
    for (uint32_t p : small){
        if ((uint64_t)p * p > hi) break;
        for (uint64_t x = max((uint64_t)p * p, ((uint64_t)lo + p - 1) / p * p); x <= hi; x += p) composite[x - lo] = 1;
    }
    vector<uint32_t> res;
    for (uint64_t x = max(lo, 2u); x <= hi; x++) if (!composite[x - lo]) res.push_back(x);
    return res;
}

int main(){
    fast_sieve();

    vector<uint32_t> small;
    vector<char> composite(1 << 16, 0);
    for (uint32_t i = 2; i < (1 << 16); i++){
        if (composite[i]) continue;
        small.push_back(i);
        for (uint32_t j = i * i; j < (1 << 16); j += i) composite[j] = 1;
    }

    auto check = [&](uint32_t lo, uint32_t hi){
        auto expected = segment(lo, hi, small);
        uint32_t* first = lower_bound(primes, primes + prime_cnt, lo);
        assert((uint32_t)(primes + prime_cnt - first) >= expected.size());
        assert(equal(expected.begin(), expected.end(), first));
        assert(first + expected.size() == primes + prime_cnt || first[expected.size()] > hi);
    };

    assert(prime_cnt == 105097565 && is_sorted(primes, primes + prime_cnt));
    check(0, 3000000);                                       /// the wheel and pre-sieved masks dominate the first blocks
    check(2147483647u - 3000000, 2147483647u);               /// the final, partial block
    for (long long it = 0; it < stress::scaled(40); it++){
        uint32_t lo = stress::rand_int(0, 2147483647LL - 300000), len = stress::rand_int(1, 300000);
        if (it % 4 == 0) lo = max(0LL, (long long)(lo / block_size * block_size) - stress::rand_int(0, 1000));  /// straddle a block boundary
        check(lo, lo + len);
    }
    return 0;
}
