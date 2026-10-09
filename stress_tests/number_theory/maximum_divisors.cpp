#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/maximum_divisors.cpp"
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

/// Every record holder has non-increasing exponents over consecutive primes, so enumerating all such products
/// and keeping each new divisor count record lists the answers, built here by full enumeration rather than pruned search
vector<pair<uint64_t, uint64_t>> highly_composite(uint64_t limit){
    const int ps[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    vector<pair<uint64_t, uint64_t>> all;
    function<void(int, int, uint64_t, uint64_t)> go = [&](int i, int max_e, uint64_t x, uint64_t cnt){
        all.push_back({x, cnt});
        for (int e = 1; i < 20 && e <= max_e && x <= limit / ps[i]; e++) go(i + 1, e, x *= ps[i], cnt * (e + 1));
    };
    go(0, 64, 1, 1);
    sort(all.begin(), all.end());

    vector<pair<uint64_t, uint64_t>> records;
    for (auto& v : all) if (records.empty() || v.second > records.back().second) records.push_back(v);
    return records;
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

    /// Past the sieve: the answer is the last record not exceeding the limit, the header table anchors the records to OEIS A002182
    auto records = highly_composite(1000000000000000000ULL);
    auto expected = [&](uint64_t limit){ return *prev(upper_bound(records.begin(), records.end(), Pair(limit, ULLONG_MAX))); };
    assert(expected(1000000) == Pair(720720, 240) && expected(1000000000000ULL) == Pair(963761198400ULL, 6720));
    assert(expected(100000000000000000ULL) == Pair(74801040398884800ULL, 64512));

    for (long long it = 0; it < stress::scaled(30); it++){
        long long limit = stress::rand_int(n, it % 3 ? 1000000000000LL : 1000000000000000000LL);
        if (it % 5 == 0) limit = expected(limit).first - stress::rand_int(0, 1);  /// on a record and just below one
        auto answer = solve(limit);
        assert(answer == expected(limit) && divisor_count(answer.first) == answer.second);
    }

    /// Highly composite numbers from the header table and OEIS A002182, one below each must return the previous one
    assert(solve(735134399) == Pair(698377680, 1280));
    assert(solve(1000000000000000000LL) == Pair(897612484786617600ULL, 103680));
    assert(solve(897612484786617599LL) == Pair(748010403988848000ULL, 98304));

    return 0;
}
