/***
 *
 * Du's Sieve
 * Prefix sums of Euler's phi and Mobius mu up to n = 1e11 in sublinear time
 *
 * Complexity: O(n^(2/3)) time with the default table of n^(2/3) / 2 entries, 12 bytes per table entry
 *
 * DuSieve du(n): sieves phi and mu up to limit = n^(2/3) / 2 and memoizes the prefix sums above it
 * DuSieve du(n, table): limit = max(1, min(n, table)), queries then cost O(n / sqrt(table)) for table >= sqrt(n)
 * du.mertens(x): mu(1) + ... + mu(x) for 0 <= x < 2^53
 * du.phi_sum(x): phi(1) + ... + phi(x) for 0 <= x < 2^53, exact, reduce it modulo a prime with % mod
 *
 * Notes:
 *   - Memory: 32 MB at n = 1e10, 134 MB at n = 1e11
 *   - phi_sum(n) plus mertens(n) at -O2: 0.14 s at n = 1e10, 0.66 s at n = 1e11 (0.22 s sieve, 0.44 s queries)
 *   - A smaller table trades time for memory, a larger one memory for time: n = 1e11, table = 2e6 takes 33 MB, 1.1 s
 *   - Queries above n cost more than O(n^(2/3))
 *   - Min25Sieve in min_25_sieve.cpp sums other multiplicative functions
 *   - Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct DuSieve{
    int limit;
    vector<long long> phi_prefix;
    vector<int> mu_prefix;
    unordered_map<long long, __int128> phi_memo;
    unordered_map<long long, long long> mu_memo;

    /// n^(2/3) / 2 entries ran 12% faster than n^(2/3) at n = 1e10 and 1e11, in half the memory
    DuSieve(long long n) : DuSieve(n, (int)min(2e9, pow((double)max(n, 1LL), 2.0 / 3.0) / 2)) {}

    DuSieve(long long n, int table) : limit(max(1, (int)min<long long>(n, table))), phi_prefix(limit + 1, 0), mu_prefix(limit + 1, 0){
        vector<int> primes;

        phi_prefix[1] = mu_prefix[1] = 1;

        /// Linear sieve writing phi and mu in place, phi_prefix[i] == 0 means no smaller i reached it, so i is prime
        for (int i = 2; i <= limit; i++){
            if (!phi_prefix[i]){
                phi_prefix[i] = i - 1, mu_prefix[i] = -1;
                primes.push_back(i);
            }
            for (int p : primes){
                long long x = (long long)i * p;
                if (x > limit) break;
                if (i % p == 0){
                    phi_prefix[x] = phi_prefix[i] * p, mu_prefix[x] = 0;
                    break;
                }
                phi_prefix[x] = phi_prefix[i] * (p - 1), mu_prefix[x] = -mu_prefix[i];
            }
        }

        for (int i = 2; i <= limit; i++){
            phi_prefix[i] += phi_prefix[i - 1];
            mu_prefix[i] += mu_prefix[i - 1];
        }
    }

    long long mertens(long long x){
        assert(0 <= x && x < (1LL << 53));
        if (x <= limit) return mu_prefix[x];
        auto it = mu_memo.find(x);
        if (it != mu_memo.end()) return it->second;

        /// Splitting at sqrt(x) and dividing in double, exact for x < 2^53, made queries 12x faster at n = 1e11
        long long r = isqrt(x), res = 1;
        double xd = x;
        for (long long i = 2; i <= r; i++){
            long long q = xd / i;
            res -= q <= limit ? mu_prefix[q] : mertens(q);
        }
        for (long long q = 1, last = x; q <= x / (r + 1); q++){
            long long next = xd / (q + 1);
            res -= (last - next) * (q <= limit ? mu_prefix[q] : mertens(q));
            last = next;
        }

        return mu_memo[x] = res;
    }

    __int128 phi_sum(long long x){
        assert(0 <= x && x < (1LL << 53));
        if (x <= limit) return phi_prefix[x];
        auto it = phi_memo.find(x);
        if (it != phi_memo.end()) return it->second;

        long long r = isqrt(x);
        double xd = x;
        __int128 res = (__int128)x * (x + 1) / 2;
        for (long long i = 2; i <= r; i++){
            long long q = xd / i;
            res -= q <= limit ? (__int128)phi_prefix[q] : phi_sum(q);
        }
        for (long long q = 1, last = x; q <= x / (r + 1); q++){
            long long next = xd / (q + 1);
            res -= (last - next) * (q <= limit ? (__int128)phi_prefix[q] : phi_sum(q));
            last = next;
        }

        return phi_memo[x] = res;
    }

    static long long isqrt(long long x){
        long long r = sqrt((double)x);
        while (r * r > x) r--;
        while ((r + 1) * (r + 1) <= x) r++;
        return r;
    }
};

int main(){
    /// Prefix sums of phi and mu, from phi(1..10) = 1 1 2 2 4 2 6 4 6 4 and mu(1..10) = 1 -1 -1 0 -1 1 -1 0 0 1
    DuSieve small(100);
    assert(small.phi_sum(0) == 0);
    assert(small.phi_sum(1) == 1);
    assert(small.phi_sum(10) == 32);
    assert(small.phi_sum(100) == 3044);
    assert(small.mertens(0) == 0);
    assert(small.mertens(1) == 1);
    assert(small.mertens(10) == -1);
    assert(small.mertens(100) == 1);
    assert(small.phi_sum(1000) == 304192);
    assert(small.mertens(1000) == 2);

    DuSieve tiny(100, 1);
    assert(tiny.limit == 1);
    assert(tiny.phi_sum(10) == 32);
    assert(tiny.phi_sum(100) == 3044);
    assert(tiny.mertens(100) == 1);

    DuSieve zero(0);
    assert(zero.phi_sum(0) == 0);
    assert(zero.mertens(0) == 0);
    assert(zero.phi_sum(2) == 2);
    assert(zero.mertens(2) == 0);

    /// Values at 10^k from OEIS A064018 and A084237, at 10^k - 1 by subtracting phi(10^k) = 0.4 * 10^k and mu(10^k) = 0
    DuSieve big(100000000000LL);
    assert(big.phi_sum(1000000000) == 303963551173008414LL);
    assert(big.mertens(1000000000) == -222);
    assert(big.phi_sum(999999999) == 303963550773008414LL);
    assert(big.mertens(100000000) == 1928);
    assert(big.mertens(10000000000LL) == -33722);
    assert(big.mertens(100000000000LL) == -87856);
    assert(big.phi_sum(10000000000LL) == (__int128)3039635509288621636LL * 10 + 6);
    assert(big.phi_sum(100000000000LL) == (__int128)3039635509283386211LL * 1000 + 140);
    assert(big.phi_sum(99999999999LL) == (__int128)3039635509243386211LL * 1000 + 140);
    assert(big.mertens(99999999999LL) == -87856);
    assert(big.phi_sum(9999999999LL) == (__int128)3039635508888621636LL * 10 + 6);
    assert(big.mertens(9999999999LL) == -33722);

    return 0;
}
