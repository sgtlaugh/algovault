/***
 *
 * Min_25 Sieve
 * Prefix sum of a multiplicative function f(1) + f(2) + ... + f(n) modulo m, from f(p) as a polynomial in p and f on prime powers
 *
 * Complexity: O(deg * n^(3/4) / log n) time, O(sqrt(n) * deg) memory, where deg = poly.size()
 * The prime sum phase meets this bound, the recursive phase is O(n^(1 - eps)) in theory and runs at the same order for n <= 1e11
 * Measured at -O2: about 0.9 s at n = 1e10 (1.2 s with a cubic f(p)), 4.7 s at n = 1e11
 *
 * Min25Sieve sieve(n, mod, poly): 1 <= n <= 1e11, 1 <= mod <= 2^32 (any modulus, prime or not)
 *     poly[k] is the coefficient of p^k in f(p) for primes p, coefficients may be negative
 * sieve.sum(f): f(1) + ... + f(n) mod m, where f(p, e, pe) returns f(p^e) for a prime p, e >= 1 and pe = p^e <= n
 *     f must agree with poly at e = 1, its result may be negative, it is reduced mod m
 *     sum can be called with several f sharing the same f(p)
 *
 * Euler phi:      Min25Sieve(n, mod, {-1, 1}).sum([](long long p, int, long long pe){ return pe / p * (p - 1); })
 * Mobius mu:      Min25Sieve(n, mod, {-1}).sum([](long long, int e, long long){ return e == 1 ? -1 : 0; })
 * Divisor count:  Min25Sieve(n, mod, {2}).sum([](long long, int e, long long){ return e + 1; })
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Min25Sieve{
    long long n;
    int r, cnt;
    uint64_t mod, inv;
    vector<int> primes;
    vector<long long> vals;                     /// every distinct n / i, descending: n / 1, n / 2, ... down to r, then r - 1, ..., 1
    vector<uint64_t> prime_sum, prime_prefix;   /// sum of f(p) over primes p <= vals[i], over the first j primes

    Min25Sieve(long long n, uint64_t mod, const vector<long long>& poly) : n(n), mod(mod){
        assert(n >= 1 && 1 <= mod && mod <= (1ULL << 32));
        inv = UINT64_MAX / mod;
        r = sqrtl((long double)n);
        while ((long long)r * r > n) r--;
        while ((long long)(r + 1) * (r + 1) <= n) r++;

        vector<char> composite(r + 1, 0);
        for (int p = 2; p <= r; p++){
            if (composite[p]) continue;
            primes.push_back(p);
            for (long long j = (long long)p * p; j <= r; j += p) composite[j] = 1;
        }

        for (long long i = 1; i <= n; i = n / (n / i) + 1) vals.push_back(n / i);
        cnt = vals.size();
        int deg = poly.size(), pc = primes.size(), large = cnt - r;

        /// surj[k][j] = j! S(k, j), the number of maps from k items onto j items
        vector<vector<uint64_t>> surj(deg, vector<uint64_t>(deg, 0));
        for (int k = 0; k < deg; k++){
            surj[k][0] = k == 0 ? 1 % mod : 0;
            for (int j = 1; j <= k; j++) surj[k][j] = mul(j % mod, add(surj[k - 1][j - 1], surj[k - 1][j]));
        }

        vector<uint64_t> coef(deg), g(cnt * deg), pref((pc + 1) * deg, 0), pw(pc * deg);
        for (int k = 0; k < deg; k++) coef[k] = reduce(poly[k]);
        for (int i = 0; i < cnt; i++) power_sums(vals[i], surj, g.data() + i * deg);
        for (int i = 0; i < pc; i++){
            uint64_t pk = 1 % mod;
            for (int k = 0; k < deg; k++){
                pw[i * deg + k] = pk;
                pref[(i + 1) * deg + k] = add(pref[i * deg + k], pk);
                pk = mul(pk, primes[i] % mod);
            }
        }

        /// g[j][k] goes from the sum of x^k over 2 <= x <= vals[j] to the sum over primes only, sifting out multiples of each p
        /// A large vals[j] is n / (j + 1), so vals[j] / p = n / d for d = (j + 1) p, which is vals[d - 1] while still large
        /// A small vals[j] is cnt - j and fits in an int, 32-bit division is several times faster than 64-bit on older judges
        for (int i = 0; i < pc; i++){
            int p = primes[i];
            for (int j = 0; j < cnt && (long long)p * p <= vals[j]; j++){
                long long d = (long long)(j + 1) * p;
                int to = j >= large ? cnt - (cnt - j) / p : d <= large ? d - 1 : cnt - n / d;
                for (int k = 0; k < deg; k++){
                    uint64_t diff = sub(g[to * deg + k], pref[i * deg + k]);
                    g[j * deg + k] = sub(g[j * deg + k], k ? mul(pw[i * deg + k], diff) : diff);
                }
            }
        }

        prime_sum.assign(cnt, 0), prime_prefix.assign(pc + 1, 0);
        for (int i = 0; i < cnt; i++){
            for (int k = 0; k < deg; k++) prime_sum[i] = add(prime_sum[i], mul(coef[k], g[i * deg + k]));
        }
        for (int i = 0; i <= pc; i++){
            for (int k = 0; k < deg; k++) prime_prefix[i] = add(prime_prefix[i], mul(coef[k], pref[i * deg + k]));
        }
    }

    template <typename F>
    uint64_t sum(F f) const{
        return add(rec(n, 0, f), 1 % mod);
    }

    uint64_t add(uint64_t a, uint64_t b) const{
        return a + b >= mod ? a + b - mod : a + b;
    }

    /// Barrett reduction of any 64-bit x, inv = floor((2^64 - 1) / m) leaves the quotient estimate at most 1 short
    uint64_t fold(uint64_t x) const{
        uint64_t res = x - (uint64_t)((unsigned __int128)x * inv >> 64) * mod;
        return res >= mod ? res - mod : res;
    }

    int index(long long v) const{
        return v <= r ? cnt - v : n / v - 1;
    }

    uint64_t mul(uint64_t a, uint64_t b) const{
        return fold(a * b);
    }

    /// out[k] = sum of x^k over 2 <= x <= v, from sum_{x=0}^{v} x^k = sum_j surj[k][j] C(v + 1, j + 1)
    /// C(v + 1, j + 1) is built exactly by dividing each factor of (j + 1)! out of the numerators, so m need not be prime
    void power_sums(long long v, const vector<vector<uint64_t>>& surj, uint64_t* out) const{
        int deg = surj.size();
        vector<long long> num;
        vector<uint64_t> binom(deg, 0);
        for (int j = 0; j < deg && j <= v; j++){
            num.push_back(v + 1 - j);
            for (int x = j + 1, q = 2; x > 1; q++){
                for (; x % q == 0; x /= q){
                    for (auto& y : num) if (y % q == 0){ y /= q; break; }
                }
            }

            binom[j] = 1 % mod;
            for (long long y : num) binom[j] = mul(binom[j], fold(y));
        }

        for (int k = 0; k < deg; k++){
            uint64_t total = 0;
            for (int j = 0; j <= k; j++) total = add(total, mul(surj[k][j], binom[j]));
            out[k] = sub(total, (1 + (k == 0)) % mod);
        }
    }

    /// sum of f(x) over 2 <= x <= v with smallest prime factor >= primes[j], requires v >= primes[j - 1]
    template <typename F>
    uint64_t rec(long long v, int j, F& f) const{
        uint64_t res = sub(prime_sum[index(v)], prime_prefix[j]);

        for (int i = j; i < (int)primes.size() && (long long)primes[i] * primes[i] <= v; i++){
            long long p = primes[i];
            int e = 1;
            for (long long pe = p; pe <= v / p; pe *= p, e++){
                res = add(res, add(mul(reduce(f(p, e, pe)), rec(v / pe, i + 1, f)), reduce(f(p, e + 1, pe * p))));
            }
        }
        return res;
    }

    uint64_t reduce(long long x) const{
        return x >= 0 ? fold(x) : sub(0, fold(-(uint64_t)x));
    }

    uint64_t sub(uint64_t a, uint64_t b) const{
        return a >= b ? a - b : a + mod - b;
    }
};

int main(){
    auto phi = [](long long p, int, long long pe){ return pe / p * (p - 1); };
    auto divisors = [](long long, int e, long long){ return (long long)e + 1; };
    assert(Min25Sieve(10, 1000, {-1, 1}).sum(phi) == 32);   /// phi(1..10) = 1 1 2 2 4 2 6 4 6 4
    assert(Min25Sieve(10, 1000, {2}).sum(divisors) == 27);  /// d(1..10) = 1 2 2 3 2 4 2 4 3 4
    assert(Min25Sieve(100, 1000000, {-1, 1}).sum(phi) == 3044);

    /// One sieve, two functions with f(p) = -1: Mobius and Liouville
    Min25Sieve minus_one(10, 1000, {-1});
    assert(minus_one.sum([](long long, int e, long long){ return e == 1 ? -1 : 0; }) == 999);  /// Mertens M(10) = -1
    assert(minus_one.sum([](long long, int e, long long){ return e % 2 ? -1 : 1; }) == 0);

    /// f(p) = p^2, so poly = {0, 0, 1}: the sum of squares
    assert(Min25Sieve(10, 1000, {0, 0, 1}).sum([](long long, int, long long pe){ return pe * pe; }) == 385);
    return 0;
}
