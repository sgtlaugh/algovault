#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/min_25_sieve.cpp"
#undef main

using Func = function<long long(long long, int, long long)>;

/// Linear sieve tables for 1..limit, built once: p^e is the exact power of the smallest prime factor p of x
struct LinearTable{
    int limit;
    vector<int> spf, power, exponent;

    LinearTable(int limit) : limit(limit), spf(limit + 1, 0), power(limit + 1, 0), exponent(limit + 1, 0){
        vector<int> primes;
        for (int i = 2; i <= limit; i++){
            if (!spf[i]) spf[i] = power[i] = i, exponent[i] = 1, primes.push_back(i);
            for (int p : primes){
                long long x = (long long)i * p;
                if (p > spf[i] || x > limit) break;
                spf[x] = p;
                power[x] = p == spf[i] ? power[i] * p : p;
                exponent[x] = p == spf[i] ? exponent[i] + 1 : 1;
            }
        }
    }

    /// f(1) + ... + f(x) mod m for every x, from f(x) = f(x / p^e) * f(p^e)
    vector<uint64_t> prefix(uint64_t mod, const Func& f) const{
        long long m = mod;
        vector<uint64_t> value(limit + 1, 0), res(limit + 1, 0);
        value[1] = 1 % mod;
        for (int x = 2; x <= limit; x++){
            if (x == power[x]) value[x] = (f(spf[x], exponent[x], x) % m + m) % m;
            else value[x] = value[x / power[x]] * value[power[x]] % mod;
        }

        for (int x = 1; x <= limit; x++) res[x] = (res[x - 1] + value[x]) % mod;
        return res;
    }
};

long long isqrt(long long n){
    long long r = sqrtl((long double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

/// Sums over all x <= n through the divisor side: sum_{x <= n} sum_{d | x} d^k = sum_{d <= n} d^k floor(n / d), in floor(n / d) blocks
uint64_t divisor_power_summatory(long long n, int k, uint64_t mod){
    auto power_prefix = [&](__int128 x) -> __int128 {
        if (k == 0) return x;
        if (k == 1) return x * (x + 1) / 2;
        return x * (x + 1) * (2 * x + 1) / 6;
    };

    __int128 total = 0;
    for (long long lo = 1, hi; lo <= n; lo = hi + 1){
        long long q = n / lo;
        hi = n / q;
        total += (power_prefix(hi) - power_prefix(lo - 1)) % mod * q % mod;
    }
    return total % mod;
}

/// Count of squarefree x <= n, sum_{d <= sqrt(n)} mu(d) floor(n / d^2)
long long squarefree_count(long long n){
    int r = isqrt(n);
    vector<int> mu(r + 1, 1);
    vector<char> composite(r + 1, 0);
    for (int p = 2; p <= r; p++){
        if (composite[p]) continue;
        for (int j = p; j <= r; j += p) composite[j] = j > p, mu[j] = -mu[j];
        for (long long j = (long long)p * p; j <= r; j += (long long)p * p) mu[j] = 0;
    }

    long long res = 0;
    for (long long d = 1; d <= r; d++) res += mu[d] * (n / (d * d));
    return res;
}

/// sum_{x <= n} x^k mod m by periodicity of x^k mod m, only for small m
uint64_t periodic_power_sum(long long n, int k, uint64_t mod){
    auto direct = [&](long long upto){
        uint64_t s = 0;
        for (long long x = 1; x <= upto; x++){
            uint64_t t = 1 % mod;
            for (int i = 0; i < k; i++) t = t * (x % mod) % mod;
            s = (s + t) % mod;
        }
        return s;
    };
    return ((uint64_t)(n / mod) % mod * direct(mod) + direct(n % mod)) % mod;
}

/// Segmented sieve up to n, checks prime_sum at every vals[i] against the running sum of f(p) over primes p <= vals[i]
void check_prime_sums(const Min25Sieve& sieve, const Func& f){
    long long n = sieve.n, m = sieve.mod;
    int i = sieve.cnt - 1;
    uint64_t total = 0;
    for (; i >= 0 && sieve.vals[i] < 2; i--) assert(sieve.prime_sum[i] == 0);
    const int block = 1 << 18;
    vector<char> composite(block);
    for (long long lo = 2; lo <= n; lo += block){
        long long hi = min(n, lo + block - 1);
        fill(composite.begin(), composite.end(), 0);
        for (long long p : sieve.primes){
            if (p * p > hi) break;
            for (long long j = max(p * p, (lo + p - 1) / p * p); j <= hi; j += p) composite[j - lo] = 1;
        }

        for (long long x = lo; x <= hi; x++){
            if (!composite[x - lo]) total = (total + (f(x, 1, x) % m + m)) % m;
            for (; i >= 0 && sieve.vals[i] == x; i--) assert(sieve.prime_sum[i] == total);
        }
    }
    assert(i == -1);
}

uint64_t random_mod(){
    static const vector<uint64_t> fixed = {1, 2, 6, 1000000, 998244353, 1000000007, 469762049, 4294967291ULL, 1ULL << 32};
    if (stress::rand_int(0, 3)) return fixed[stress::rand_int(0, fixed.size() - 1)];
    return stress::rand_int(1, 1LL << 32);
}

/// A multiplicative function with f(p) = poly(p) and pseudo-random values on higher prime powers, poly has min_size to 6 terms
pair<vector<long long>, Func> random_function(uint64_t mod, int min_size = 0){
    vector<long long> poly(stress::rand_int(min_size, 6));
    for (auto& c : poly) c = stress::rand_int(-1000000000000LL, 1000000000000LL);
    long long salt = stress::rand_int(0, 1LL << 40), m = mod;

    Func f = [poly, salt, m](long long p, int e, long long pe) -> long long {
        if (e > 1) return (pe * 2654435761LL + e * 40503 + salt) % 1000000007 - 500000000;
        __int128 v = 0;
        for (int k = (int)poly.size() - 1; k >= 0; k--) v = (v * p + poly[k]) % m;
        return (long long)v;
    };
    return {poly, f};
}

int main(){
    const int limit = 1000000;
    Func phi = [](long long p, int, long long pe){ return pe / p * (p - 1); };
    Func mu = [](long long, int e, long long){ return e == 1 ? -1LL : 0LL; };
    Func divisors = [](long long, int e, long long){ return (long long)e + 1; };
    Func sigma2 = [](long long p, int, long long pe){ return (long long)(((__int128)pe * pe * p * p - 1) / (p * p - 1) % 1000000007); };
    vector<pair<vector<long long>, Func>> named = {{{-1, 1}, phi}, {{-1}, mu}, {{2}, divisors}, {{1, 0, 1}, sigma2}};

    vector<long long> queries;
    for (long long n = 1; n <= 300; n++) queries.push_back(n);
    for (long long p : {2, 3, 7, 31, 97, 997}){
        for (long long d = -1; d <= 1; d++) queries.push_back(p * p + d);
    }
    for (long long n : {65535, 65536, 999999, 1000000}) queries.push_back(n);

    /// Barrett reduction against %, over the full 64-bit range including its top
    for (long long it = 0; it < stress::scaled(300); it++){
        uint64_t mod = it < 2 ? it + (1ULL << 32) - 1 : random_mod();
        Min25Sieve sieve(1, mod, {});
        for (int t = 0; t < 3000; t++){
            uint64_t x = t < 3 ? UINT64_MAX - t : stress::rng()() >> stress::rand_int(0, 63);
            assert(sieve.fold(x) == x % mod);
        }
    }

    /// Linear sieve prefix sums up to 1e6: every n <= 300, around prime squares and random n, named and random f of degree up to 5
    LinearTable table(limit);
    for (long long it = 0; it < stress::scaled(30); it++){
        uint64_t mod = random_mod();
        auto [poly, f] = it < (long long)named.size() * 3 ? named[it % named.size()] : random_function(mod);
        if (it < (long long)named.size() * 3 && it % named.size() == 3) mod = 1000000007;

        auto prefix = table.prefix(mod, f);
        vector<long long> ns = queries;
        for (int t = 0; t < 30; t++) ns.push_back(stress::rand_int(1, limit));
        for (long long n : ns) assert(Min25Sieve(n, mod, poly).sum(f) == prefix[n]);
    }

    /// Large n against O(sqrt n) formulas: divisor count, sigma, sigma_2 and the squarefree count, plus the n = 1e11 boundary
    Func sigma = [](long long p, int, long long pe){ return (pe * p - 1) / (p - 1); };
    Func squarefree = [](long long, int e, long long){ return (long long)(e == 1); };
    for (long long it = 0; it < stress::scaled(4); it++){
        long long n = stress::rand_int(100000000LL, 1000000000LL);
        uint64_t mod = random_mod();
        assert(Min25Sieve(n, mod, {2}).sum(divisors) == divisor_power_summatory(n, 0, mod));
        assert(Min25Sieve(n, mod, {1, 1}).sum(sigma) == divisor_power_summatory(n, 1, mod));
        assert(Min25Sieve(n, 1000000007, {1, 0, 1}).sum(sigma2) == divisor_power_summatory(n, 2, 1000000007));
        assert(Min25Sieve(n, 1ULL << 32, {1}).sum(squarefree) == (uint64_t)squarefree_count(n) % (1ULL << 32));
    }

    assert(Min25Sieve(10000000000LL, 1000000007, {1, 0, 1}).sum(sigma2) == divisor_power_summatory(10000000000LL, 2, 1000000007));
    assert(Min25Sieve(100000000000LL, 1ULL << 32, {2}).sum(divisors) == divisor_power_summatory(100000000000LL, 0, 1ULL << 32));

    /// Completely multiplicative x^k up to deg 5 under composite moduli, where power sums cannot divide by k + 1
    for (uint64_t mod : {1000000ULL, 1ULL << 20, 720720ULL}){
        for (int k = 0; k <= 5; k++){
            long long n = stress::rand_int(10000000LL, 100000000LL);
            vector<long long> poly(k + 1, 0);
            poly[k] = 1;
            Func power = [k, mod](long long, int, long long pe){
                uint64_t t = 1 % mod;
                for (int i = 0; i < k; i++) t = t * (pe % mod) % mod;
                return (long long)t;
            };
            assert(Min25Sieve(n, mod, poly).sum(power) == periodic_power_sum(n, k, mod));
        }
    }

    /// Random f(p) of degree 3 to 5 with mixed coefficients at large n, the prime sum at every vals[i] against a segmented sieve
    /// The modulus stays above 2^31, a tiny one such as 1 or 2 would zero out the high coefficients this block is here to test
    for (long long it = 0; it < stress::scaled(1); it++){
        uint64_t mod = stress::rand_int(1LL << 31, 1LL << 32);
        auto [poly, f] = random_function(mod, 4);
        check_prime_sums(Min25Sieve(stress::rand_int(100000000LL, 200000000LL), mod, poly), f);
    }

    return 0;
}
