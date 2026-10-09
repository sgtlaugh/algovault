/***
 *
 * Du's Sieve
 * Prefix sums of Euler's phi and Mobius mu up to n = 1e11 in sublinear time
 *
 * Complexity: O(n^(2/3)) time with the default table of about n^(2/3) entries, 12 bytes per table entry
 * n = 1e10: 56 MB, n = 1e11: 260 MB. Requires __int128
 * Measured at -O2 for phi_sum(n) plus mertens(n): n = 1e10 takes 1.0 s, n = 1e11 takes 4.7 s (0.5 s sieve, 4.2 s queries)
 *
 * DuSieve du(n): sieves phi and mu up to limit = n^(2/3) and memoizes the prefix sums above it
 * DuSieve du(n, table): limit = max(1, min(n, table)) instead, queries then cost O(n / sqrt(table)) for table >= sqrt(n),
 *     so a smaller table trades time for memory and a larger one memory for time (n = 1e11, table = 5e7: 600 MB, 3.7 s)
 * du.phi_sum(x): phi(1) + ... + phi(x) for x >= 0, exact, reduce it modulo a prime with % mod
 * du.mertens(x): mu(1) + ... + mu(x) for x >= 0
 * Queries x <= n keep the complexity, larger x work too but cost more than O(n^(2/3))
 *
 * Both come from the Dirichlet convolutions phi * 1 = id and mu * 1 = e:
 *     sum_{d=1..x} phi_sum(x / d) = x (x + 1) / 2 and sum_{d=1..x} mertens(x / d) = 1
 * so each value above the table is its convolution total minus O(sqrt(x)) blocks of equal x / d
 *
 * Library Checker sum_of_totient_function, n <= 1e10:
 *     DuSieve du(n);
 *     long long answer = du.phi_sum(n) % 998244353;
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

    DuSieve(long long n) : DuSieve(n, (int)min(2e9, pow((double)max(n, 1LL), 2.0 / 3.0))) {}

    DuSieve(long long n, int table) : limit(max(1, (int)min<long long>(n, table))), phi_prefix(limit + 1, 0), mu_prefix(limit + 1, 0){
        vector<int> primes;
        phi_prefix[1] = mu_prefix[1] = 1;

        /// Linear sieve writing phi and mu in place before the prefix sums, phi_prefix[i] == 0 means no smaller i reached it, so i is prime
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
        assert(x >= 0);
        if (x <= limit) return mu_prefix[x];
        auto it = mu_memo.find(x);
        if (it != mu_memo.end()) return it->second;

        long long res = 1;
        for (long long i = 2, last; i <= x; i = last + 1){
            long long q = x / i;
            last = x / q;
            res -= (last - i + 1) * (q <= limit ? mu_prefix[q] : mertens(q));
        }

        return mu_memo[x] = res;
    }

    __int128 phi_sum(long long x){
        assert(x >= 0);
        if (x <= limit) return phi_prefix[x];
        auto it = phi_memo.find(x);
        if (it != phi_memo.end()) return it->second;

        __int128 res = (__int128)x * (x + 1) / 2;
        for (long long i = 2, last; i <= x; i = last + 1){
            long long q = x / i;
            last = x / q;
            res -= (last - i + 1) * (q <= limit ? (__int128)phi_prefix[q] : phi_sum(q));
        }

        return phi_memo[x] = res;
    }
};

int main(){
    DuSieve small(100);
    assert(small.phi_sum(0) == 0 && small.phi_sum(1) == 1 && small.phi_sum(10) == 32 && small.phi_sum(100) == 3044);
    assert(small.mertens(0) == 0 && small.mertens(1) == 1 && small.mertens(10) == -1 && small.mertens(100) == 1);
    assert(small.phi_sum(1000) == 304192 && small.mertens(1000) == 2);

    DuSieve tiny(100, 1);
    assert(tiny.limit == 1 && tiny.phi_sum(10) == 32 && tiny.phi_sum(100) == 3044 && tiny.mertens(100) == 1);

    DuSieve zero(0);
    assert(zero.phi_sum(0) == 0 && zero.mertens(0) == 0 && zero.phi_sum(2) == 2 && zero.mertens(2) == 0);

    /// Values at 10^k from OEIS A064018 and A084237, at 10^k - 1 by subtracting phi(10^k) = 0.4 * 10^k and mu(10^k) = 0
    DuSieve big(100000000000LL);
    assert(big.phi_sum(1000000000) == 303963551173008414LL && big.mertens(1000000000) == -222);
    assert(big.phi_sum(999999999) == 303963550773008414LL && big.mertens(100000000) == 1928);
    assert(big.mertens(10000000000LL) == -33722 && big.mertens(100000000000LL) == -87856);
    assert(big.phi_sum(10000000000LL) == (__int128)3039635509288621636LL * 10 + 6);
    assert(big.phi_sum(100000000000LL) == (__int128)3039635509283386211LL * 1000 + 140);
    assert(big.phi_sum(99999999999LL) == (__int128)3039635509243386211LL * 1000 + 140 && big.mertens(99999999999LL) == -87856);
    assert(big.phi_sum(9999999999LL) == (__int128)3039635508888621636LL * 10 + 6 && big.mertens(9999999999LL) == -33722);

    return 0;
}
