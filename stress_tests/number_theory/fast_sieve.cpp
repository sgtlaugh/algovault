#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/fast_sieve.cpp"
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
    vector<uint32_t> small;
    vector<char> composite(1 << 16, 0);
    for (uint32_t i = 2; i < (1 << 16); i++){
        if (composite[i]) continue;
        small.push_back(i);
        for (uint32_t j = i * i; j < (1 << 16); j += i) composite[j] = 1;
    }

    /// A small n must stop exactly at n, including n inside a word, at a block boundary and in a later block
    const uint32_t block_size = FastSieve::block_size;
    for (uint32_t n = 0; n <= 40; n++) assert(FastSieve(n).primes == segment(0, n, small));
    for (long long it = 0; it < stress::scaled(20); it++){
        uint32_t n = stress::rand_int(41, 300000);
        assert(FastSieve(n).primes == segment(0, n, small));
    }
    for (uint32_t n : {block_size - 1, block_size, block_size + 1, 3 * block_size + 12345, 5 * block_size - 2}){
        assert(FastSieve(n).primes == segment(0, n, small));
    }
    for (long long it = 0; it < stress::scaled(3); it++){
        uint32_t n = stress::rand_int(1, 20000000);
        assert(FastSieve(n).primes == segment(0, n, small));
    }

    FastSieve sieve(2147483647);
    const vector<uint32_t>& primes = sieve.primes;

    auto check = [&](uint32_t lo, uint32_t hi){
        auto expected = segment(lo, hi, small);
        auto first = lower_bound(primes.begin(), primes.end(), lo);
        assert((uint32_t)(primes.end() - first) >= expected.size());
        assert(equal(expected.begin(), expected.end(), first));
        assert(first + expected.size() == primes.end() || first[expected.size()] > hi);
    };

    assert(primes.size() == 105097565 && is_sorted(primes.begin(), primes.end()));
    check(0, 3000000);                                       /// the wheel and pre-sieved masks dominate the first blocks
    check(2147483647u - 3000000, 2147483647u);               /// the final, partial block
    for (long long it = 0; it < stress::scaled(40); it++){
        uint32_t lo = stress::rand_int(0, 2147483647LL - 300000), len = stress::rand_int(1, 300000);
        if (it % 4 == 0) lo = max(0LL, (long long)(lo / block_size * block_size) - stress::rand_int(0, 1000));  /// straddle a block boundary
        check(lo, lo + len);
    }

    return 0;
}
