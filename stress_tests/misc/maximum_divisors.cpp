#include "../common.h"

#define main library_main
#include "../../code_library/misc/maximum_divisors.cpp"
#undef main

/// Any candidate is a product of the first primes, so trial division by those alone must reduce it to 1
uint64_t divisor_count(uint64_t x){
    uint64_t res = 1;
    for (uint64_t p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71}){
        int e = 0;
        while (x % p == 0) x /= p, e++;
        res *= e + 1;
    }
    assert(x == 1);
    return res;
}

int main(){
    typedef pair<uint64_t, uint64_t> Pair;
    const int n = 2000000;
    vector<int> d(n + 1, 0);
    for (int i = 1; i <= n; i++) for (int j = i; j <= n; j += i) d[j]++;

    /// best[x]: the smallest number up to x with the most divisors
    vector<int> best(n + 1, 1);
    for (int x = 2; x <= n; x++) best[x] = d[x] > d[best[x - 1]] ? x : best[x - 1];

    for (int x = 1; x <= 3000; x++) assert(solve(x) == Pair(best[x], d[best[x]]));
    for (long long it = 0; it < stress::scaled(3000); it++){
        int x = stress::rand_int(1, n);
        assert(solve(x) == Pair(best[x], d[best[x]]));
    }

    /// Past the sieve: the answer must have the claimed count, and nothing smaller may match it, so answer - 1 must fall short
    for (long long it = 0; it < stress::scaled(30); it++){
        long long limit = stress::rand_int(n, it % 3 ? 1000000000000LL : 1000000000000000000LL);
        auto [x, cnt] = solve(limit);
        assert(x <= (uint64_t)limit && divisor_count(x) == cnt);
        assert(solve(x) == Pair(x, cnt) && solve(x - 1).second < cnt);
    }

    /// Highly composite numbers from the header table and OEIS A002182, one below each must return the previous one
    assert(solve(735134399) == Pair(698377680, 1280));
    assert(solve(1000000000000000000LL) == Pair(897612484786617600ULL, 103680));
    assert(solve(897612484786617599LL) == Pair(748010403988848000ULL, 98304));
    return 0;
}
