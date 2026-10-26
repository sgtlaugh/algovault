/***
 *
 * Fast Prime Sums
 * Sum of primes not exceeding n in sublinear time with the Meissel-Lehmer algorithm
 *
 * Complexity: O(maxv + maxn maxm) to build, a query has no tight bound and grows about 10x per decade of n
 * Memory: 12 * maxv + 8 * maxn * maxm bytes, about 95 MB with the defaults (n < 2.5e13)
 *
 * PrimeSum ps(maxv, maxn, maxm): sieves below maxv and tabulates phi for the first maxn primes and m < maxm
 *     maxv must exceed the maxn-th prime
 * ps.sum(n): the sum of primes not exceeding n, for 0 <= n < maxv^2
 *
 * Notes:
 *   - A query with the defaults at -O2: 0.03 s at n = 1e11, 0.23 s at 1e12, 2.2 s at 1e13
 *   - PrimeSum ps(20000010, 50, 2000010) takes about 1 GB and 0.8 s to build, then queries near 1e13 run 3x faster
 *   - Requires __int128
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// maxn and maxm stay runtime arguments: as template parameters, like PrimeCounter's, sum(1e13) took 2.07 s vs 2.16 s
struct PrimeSum{
    int maxv, maxn, maxm;
    vector<int> primes, pi;
    vector<uint64_t> pi_sum, dp;   /// dp[i * maxm + j] = phi(j, i)

    PrimeSum(int maxv = 5000000, int maxn = 20, int maxm = 200000) : maxv(maxv), maxn(maxn), maxm(maxm){
        sieve();
        assert(maxn >= 1 && maxm >= 1 && (int)primes.size() >= maxn);

        dp.resize((size_t)maxn * maxm);
        for (int j = 0; j < maxm; j++) dp[j] = (uint64_t)j * (j + 1) / 2;
        for (int i = 1; i < maxn; i++){
            uint64_t *cur = &dp[(size_t)i * maxm], *prev = cur - maxm;
            for (int j = 1; j < maxm; j++){
                cur[j] = prev[j] - prev[fast_div(j, primes[i - 1])] * primes[i - 1];
            }
        }
    }

    /// uint64_t only below 2^32, where m (m + 1) in phi(m, 0) fits before halving and the rest wraps exactly mod 2^64
    __int128 sum(long long n){
        assert(0 <= n && n < (long long)maxv * maxv);
        if (n <= UINT_MAX) return lehmer((uint64_t)n);
        return lehmer((__int128)n);
    }

    void sieve(){
        vector<bool> is_prime(maxv);

        is_prime[2] = true;
        for (int i = 3; i < maxv; i += 2) is_prime[i] = true;

        for (int i = 3; i * i < maxv; i += 2){
            if (!is_prime[i]) continue;
            for (int j = i * i; j < maxv; j += (i << 1)) is_prime[j] = false;
        }

        pi.resize(maxv), pi_sum.resize(maxv);
        int cnt = 0;
        uint64_t total = 0;
        for (int i = 1; i < maxv; i++){
            if (is_prime[i]) cnt++, total += i, primes.push_back(i);
            pi[i] = cnt, pi_sum[i] = total;
        }
    }

    static uint64_t fast_div(uint64_t a, uint32_t b){
        return double(a) / b + 1e-9;
    }

    /// floor(n^(1/k)) for k = 2 or 3: the double root alone is one too low at some cubes above 1.5e15, like 114502^3
    static int root(long long n, int k){
        long long r = k == 2 ? sqrt(0.5 + n) : cbrt(0.5 + n);
        auto power = [k](long long x){ return k == 2 ? x * x : x * x * x; };
        while (power(r) > n) r--;
        while (power(r + 1) <= n) r++;
        return r;
    }

    /// Sub-sums go through sum(), so a quotient of 2^32 or more never runs in uint64_t
    template <typename T>
    T lehmer(T n){
        if (n < (T)maxv) return pi_sum[n];

        int s = root(n, 2), c = root(n, 3);
        T res = phi(n, pi[c]) + pi_sum[c] - 1;
        for (int i = pi[c]; i < pi[s]; i++){
            long long q = fast_div(n, primes[i]);
            /// pi_sum[q] read here rather than in sum(): 2.16 s vs 2.36 s for sum(1e13) at -O2
            res -= (T)((q < maxv ? pi_sum[q] : sum(q)) - pi_sum[primes[i] - 1]) * primes[i];
        }

        return res;
    }

    template <typename T>
    T phi(T m, int n){
        if (!n) return m * (m + 1) / 2;
        if (n < maxn && m < (T)maxm) return dp[(size_t)n * maxm + (size_t)m];

        uint64_t p = primes[n - 1];
        if (m < (T)maxv && p * p >= m) return pi_sum[m] - pi_sum[p] + 1;
        if (m < (T)maxv && p * p * p >= m){
            int lim = pi[(int)sqrt(0.5 + m)];
            T res = pi_sum[m] - pi_sum[p] + 1;
            for (int i = n; i < lim; i++) res += (T)primes[i] * (pi_sum[fast_div(m, primes[i])] - pi_sum[primes[i] - 1]);
            return res;
        }

        return phi(m, n - 1) - phi((T)fast_div(m, p), n - 1) * p;
    }
};

int main(){
    /// The default tables answer n < 2.5e13
    PrimeSum ps;
    assert(ps.sum(0) == 0);
    assert(ps.sum(10) == 17);
    assert(ps.sum(1000) == 76127);
    assert(ps.sum(1000000) == 37550402023LL);
    assert(ps.sum(1000000000000LL) == (__int128)15929208151LL * 1157344946327LL);  /// 18435588552550705911377

    /// Small tables answer n up to maxv^2 - 1
    PrimeSum small(1000, 5, 100);
    assert(small.sum(999999) == 37550402023LL);

    return 0;
}
