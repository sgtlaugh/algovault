/***
 *
 * Prime sum function in sublinear time with the Meissel-Lehmer algorithm
 *
 * PrimeSum ps(maxv, maxn, maxm): sieves below maxv and tabulates phi for the first maxn primes and m < maxm
 *     maxv must exceed the maxn-th prime, the defaults are tuned for n around 1e13
 * ps.prime_sum(n): the sum of primes not exceeding n, for 0 <= n < maxv^2
 * It is just a templatized wrapper of lehmer(n)
 *
 * Complexity: Roughly ~O(n^(2/3)) per query, O(maxv + maxn * maxm) to build
 * Memory: 12 * maxv + 8 * maxn * maxm bytes, about 1 GB with the defaults
 * For 256 MB judges use PrimeSum ps(5000000, 20, 200000): 95 MB, about 1.9x slower on the self-test
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct PrimeSum{
    int maxv, maxn, maxm;
    vector<int> pi, primes;
    vector<uint64_t> pi_sum, dp;   /// dp[i * maxm + j] = phi(j, i)

    PrimeSum(int maxv = 20000010, int maxn = 50, int maxm = 2000010) : maxv(maxv), maxn(maxn), maxm(maxm){
        sieve();
        assert(maxn >= 1 && maxm >= 1 && (int)primes.size() > maxn);

        dp.resize((size_t)maxn * maxm);
        for (int j = 0; j < maxm; j++) dp[j] = (uint64_t)j * (j + 1) / 2;
        for (int i = 1; i < maxn; i++){
            uint64_t *cur = &dp[(size_t)i * maxm], *prev = cur - maxm;
            for (int j = 1; j < maxm; j++){
                cur[j] = prev[j] - prev[fast_div(j, primes[i])] * primes[i];
            }
        }
    }

    __int128 prime_sum(long long n){
        if (n <= UINT_MAX) return lehmer((uint64_t)n);
        return lehmer((__int128)n);
    }

private:
    static uint64_t fast_div(uint64_t a, uint32_t b){
        return double(a) / b + 1e-9;
    }

    template <typename T>
    T phi(T m, int n){
        if (!n) return (T)m * (m + 1) / 2;
        if (n < maxn && m < (T)maxm) return dp[(size_t)n * maxm + (size_t)m];
        if (m < (T)maxv && (uint64_t)primes[n] * primes[n] >= m) return pi_sum[m] - pi_sum[primes[n]] + 1;
        return phi(m, n - 1) - phi((T)fast_div(m, primes[n]), n - 1) * primes[n];
    }

    template <typename T>
    T lehmer(T n){
        if (n < (T)maxv) return pi_sum[n];

        int s = sqrt(0.5 + n), c = cbrt(0.5 + n);
        T res = phi(n, pi[c]) + pi_sum[c] - 1;

        for (int i = pi[c] + 1; i <= pi[s]; i++){
            T w = lehmer(fast_div(n, primes[i])) - pi_sum[primes[i] - 1];
            res -= w * primes[i];
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

        primes.push_back(-1);
        pi.resize(maxv), pi_sum.resize(maxv);
        int cnt = 0;
        uint64_t sum = 0;
        for (int i = 1; i < maxv; i++){
            if (is_prime[i]) cnt++, sum += i, primes.push_back(i);
            pi[i] = cnt, pi_sum[i] = sum;
        }
    }
};

int main(){
    auto start = clock();
    PrimeSum ps;
    fprintf(stderr, "Pre-process time = %0.3f\n\n", (clock()-start) / (double)CLOCKS_PER_SEC);  /// 0.940

    start = clock();

    assert(ps.prime_sum(1000) == 76127);
    assert(ps.prime_sum(1000000) == 37550402023);
    assert(ps.prime_sum(1e9) == 24739512092254535LL);
    assert(ps.prime_sum(1e10) == 2220822432581729238LL);
    assert(ps.prime_sum(UINT_MAX) == 425649736193687430LL);

    assert(ps.prime_sum(1e11) == (__int128)603698 * 333721625289043LL);      /// 201467077743744681014
    assert(ps.prime_sum(1e12) == (__int128)15929208151LL * 1157344946327LL); /// 18435588552550705911377
    assert(ps.prime_sum(1e13) == (__int128)10166702 * 167138413556114797LL); /// 1699246443377779418889494

    fprintf(stderr, "\nCalculation time = %0.3f\n", (clock()-start) / (double)CLOCKS_PER_SEC);  /// 2.130
    return 0;
}
