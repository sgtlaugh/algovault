/***
 *
 * Fast Prime Counting
 * Prime counting function in sublinear time with the Meissel-Lehmer algorithm
 *
 * Complexity: O(maxv + MAXN MAXM) to build, a query has no tight bound and grows about 10x per decade of n
 * Memory: 4 * (maxv + MAXN * MAXM) bytes, about 190 MB with the defaults
 *
 * PrimeCounter<MAXN, MAXP> pc(maxv): sieves below maxv and tabulates phi for the first MAXN primes
 *     MAXM, the phi table width, is the product of the first MAXP primes, 1 <= MAXP <= min(MAXN, 8)
 *     maxv must exceed the MAXN-th prime, the defaults are tuned for n around 1e12
 * pc.count(n): the number of primes not exceeding n, for 0 <= n < maxv^2
 *
 * Notes:
 *   - A query with the defaults at -O2: 0.006 s at n = 1e11, 0.05 s at 1e12, 0.44 s at 1e13, 4.3 s at 1e14
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// MAXN and MAXP are template parameters: as runtime members count(1e13) took 0.508 s vs 0.480 s at -O2
template<int MAXN = 50, int MAXP = 7>
struct PrimeCounter{
    static_assert(1 <= MAXP && MAXP <= MAXN && MAXP <= 8);
    static constexpr int MAXM = [](){
        int res = 1, small[] = {2, 3, 5, 7, 11, 13, 17, 19};
        for (int i = 0; i < MAXP; i++) res *= small[i];
        return res;
    }();

    int maxv;
    vector<int> primes, prod, pi, dp;   /// dp[i * MAXM + j] = phi(j, i)

    PrimeCounter(int maxv = 20000010) : maxv(maxv){
        sieve();
        assert((int)primes.size() >= MAXN);

        prod.resize(MAXP);
        prod[0] = primes[0];
        for (int i = 1; i < MAXP; i++) prod[i] = prod[i - 1] * primes[i];

        dp.resize((size_t)MAXN * MAXM);
        for (int j = 0; j < MAXM; j++) dp[j] = j;
        for (int i = 1; i < MAXN; i++){
            int *cur = &dp[(size_t)i * MAXM], *prev = cur - MAXM;
            for (int j = 1; j < MAXM; j++){
                cur[j] = prev[j] - prev[fast_div(j, primes[i - 1])];
            }
        }
    }

    uint64_t count(long long n){
        assert(0 <= n && n < (long long)maxv * maxv);
        if (n < maxv) return pi[n];

        int s = root(n, 2), c = root(n, 3);
        uint64_t res = phi(n, pi[c]) + pi[c] - 1;
        for (int i = pi[c]; i < pi[s]; i++){
            long long q = fast_div(n, primes[i]);
            /// pi[q] read here, since gcc does not inline the same check in count(): 0.45 s vs 0.52 s for count(1e13)
            res -= (q < maxv ? pi[q] : count(q)) - i;
        }

        return res;
    }

    void sieve(){
        vector<bool> is_prime(maxv);

        is_prime[2] = true;
        for (int i = 3; i < maxv; i += 2) is_prime[i] = true;

        for (int i = 3; i * i < maxv; i += 2){
            if (!is_prime[i]) continue;
            for (int j = i * i; j < maxv; j += (i << 1)) is_prime[j] = false;
        }

        pi.resize(maxv);
        for (int i = 1, cnt = 0; i < maxv; i++){
            if (is_prime[i]) cnt++, primes.push_back(i);
            pi[i] = cnt;
        }
    }

    static long long fast_div(long long a, int b){
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

    uint64_t phi(long long m, int n){
        if (!n) return m;
        if (n < MAXN && m < MAXM) return dp[(size_t)n * MAXM + m];
        if (n < MAXP) return dp[(size_t)n * MAXM + m % prod[n - 1]] + fast_div(m, prod[n - 1]) * dp[(size_t)n * MAXM + prod[n - 1]];

        long long p = primes[n - 1];
        if (m < maxv && p * p >= m) return pi[m] - n + 1;
        if (p * p * p < m || m >= maxv) return phi(m, n - 1) - phi(fast_div(m, p), n - 1);

        int lim = pi[(int)sqrt(0.5 + m)];
        uint64_t res = pi[m] - (lim + n - 2) * (lim - n + 1) / 2;
        for (int i = n; i < lim; i++){
            res += pi[fast_div(m, primes[i])];
        }

        return res;
    }
};

int main(){
    /// The default tables, tuned for n around 1e12
    PrimeCounter<> pc;
    assert(pc.count(0) == 0);
    assert(pc.count(7) == 4);
    assert(pc.count(1000) == 168);
    assert(pc.count(1000000000000LL) == 37607912018LL);

    /// Small tables answer n up to maxv^2 - 1
    PrimeCounter<10, 3> small(1000);
    assert(small.count(1000) == 168);
    assert(small.count(999999) == 78498);

    return 0;
}
