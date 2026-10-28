/***
 *
 * Sum of Two Squares
 * Every way to write n as x^2 + y^2 with x, y >= 0, for 0 <= n < 2^63
 *
 * Complexity: expected O(n^(1/4)) per Pollard rho split, O(log p) per prime p dividing n, O(k log k) for k pairs
 *
 * TwoSquares squares: holds a Pollard rho factorizer, microseconds to construct
 * squares.pairs(n): every (x, y) with x, y >= 0 and x^2 + y^2 = n, ascending, both (x, y) and (y, x) when x != y
 * squares.r2(n): number of integer pairs (x, y) with x^2 + y^2 = n, signs and order counted, 1 for n = 0
 *
 * Notes:
 *   - n has a pair exactly when every prime 3 mod 4 divides it an even number of times
 *   - At most 4608 pairs below 2^63
 *   - Requires __int128 (64-bit GCC or Clang)
 *   - Embeds checked copies of Montgomery and prm::is_prime from miller_rabin.cpp and PollardRho from pollard_rho.cpp
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Montgomery form for a fixed odd modulus n, values are stored as x * 2^64 mod n, 8x faster than __int128 % n here
/// BEGIN COPY montgomery from code_library/number_theory/miller_rabin.cpp
struct Montgomery{
    unsigned long long n, inv, r2;

    Montgomery(unsigned long long n) : n(n), inv(1){
        for (int i = 0; i < 6; i++) inv *= 2 - n * inv;  /// each Newton step doubles the correct bits of n^-1 mod 2^64
        r2 = -n % n;
        r2 = (unsigned __int128)r2 * r2 % n;
    }

    unsigned long long reduce(unsigned __int128 x) const{
        unsigned long long q = (unsigned long long)x * inv, m = ((unsigned __int128)q * n) >> 64, h = x >> 64;
        return h >= m ? h - m : h + n - m;
    }

    unsigned long long mul(unsigned long long x, unsigned long long y) const{
        return reduce((unsigned __int128)x * y);
    }

    unsigned long long add(unsigned long long x, unsigned long long y) const{
        return (x += y) >= n ? x - n : x;
    }

    unsigned long long to_mont(unsigned long long x) const{
        return mul(x, r2);
    }

    unsigned long long pow(unsigned long long x, unsigned long long e) const{
        unsigned long long res = to_mont(1);
        for (; e; e >>= 1, x = mul(x, x)){
            if (e & 1) res = mul(res, x);
        }
        return res;
    }
};
/// END COPY montgomery

/// BEGIN COPY miller_rabin from code_library/number_theory/miller_rabin.cpp
namespace prm{
    /// Jaeschke 1993, deterministic below 4759123141
    const array<int, 3> BASES_32 = {2, 7, 61};
    /// Jim Sinclair 2011, deterministic below 2^64: https://miller-rabin.appspot.com/
    const array<int, 7> BASES_64 = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

    struct SmallPrime{ unsigned long long p, inv, max_quotient; };

    /// Holds the bases 7 and 61, so n never divides a base and no base becomes 0 in Montgomery form
    const vector<SmallPrime> SMALL_PRIMES = [](){
        vector<SmallPrime> res;
        for (unsigned long long p: {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79}){
            res.push_back({p, Montgomery(p).inv, ~0ULL / p});
        }
        return res;
    }();

    /// x = a^d for odd d with n - 1 = d 2^s: a proves n composite unless x is 1 or reaches -1 within s - 1 squarings
    bool is_witness(unsigned long long x, int s, const Montgomery& mont){
        unsigned long long one = mont.to_mont(1), minus_one = mont.n - one;

        if (x == one) return false;

        for (int i = 0; i < s; i++){
            if (x == minus_one) return false;
            x = mont.mul(x, x);
        }

        return true;
    }

    template<size_t K>
    bool miller_rabin(long long n, const array<int, K>& bases){
        if (n < 2 || (n & 1) == 0) return n == 2;

        /// p | n exactly when n * p^-1 mod 2^64 <= (2^64 - 1) / p, no division: 27% faster on random n below 2^63
        for (auto& [p, inv, max_quotient]: SMALL_PRIMES){
            if ((unsigned long long)n * inv <= max_quotient) return (unsigned long long)n == p;
        }

        int s = __builtin_ctzll(n - 1);
        unsigned long long d = (n - 1) >> s;
        Montgomery mont(n);

        /// bases[0] is 2 in both lists, which alone rejects almost every composite, so the rest mostly run on primes
        if (is_witness(mont.pow(mont.to_mont(bases[0]), d), s, mont)) return false;

        /// The other bases in lockstep so their independent multiplications overlap: 63-bit primes 29% faster
        const int k = K - 1;  /// constant, so the lanes unroll: 31-bit primes 12% faster than a runtime count
        unsigned long long x[k], res[k];
        for (int i = 0; i < k; i++) x[i] = mont.to_mont(bases[i + 1]), res[i] = mont.to_mont(1);

        for (; d; d >>= 1){
            if (d & 1){
                for (int i = 0; i < k; i++) res[i] = mont.mul(res[i], x[i]);
            }
            for (int i = 0; i < k; i++) x[i] = mont.mul(x[i], x[i]);
        }

        for (int i = 0; i < k; i++){
            if (is_witness(res[i], s, mont)) return false;
        }

        return true;
    }

    bool is_prime(long long n){
        if (n < 4759123141LL) return miller_rabin(n, BASES_32);
        return miller_rabin(n, BASES_64);
    }
}
/// END COPY miller_rabin

/// BEGIN COPY pollard_rho from code_library/number_theory/pollard_rho.cpp
struct PollardRho{
    static constexpr int TRIAL_LIMIT = 1000;  /// not a 10^6 sieve: 3-7 ms less to build, 10^6 n below 10^6 25% slower
    vector<prm::SmallPrime> odd_primes;
    mt19937_64 rng;

    PollardRho() : rng(20260105){
        for (unsigned long long p = 3; p < TRIAL_LIMIT; p += 2){
            bool prime = true;
            for (unsigned long long q = 3; q * q <= p; q += 2) prime &= p % q != 0;
            if (prime) odd_primes.push_back({p, Montgomery(p).inv, ~0ULL / p});
        }
    }

    vector<long long> factorize(long long n){
        assert(n >= 1);
        vector<long long> res;
        for (; !(n & 1); n >>= 1) res.push_back(2);

        /// p | n exactly when n * p^-1 mod 2^64 <= (2^64 - 1) / p, and then n * p^-1 is n / p
        for (auto& [p, inv, max_quotient]: odd_primes){
            if ((unsigned long long)n < p * p) break;  /// n is 1 or a prime: 10^6 numbers below 10^6 factor 2x faster
            for (; (unsigned long long)n * inv <= max_quotient; n *= inv) res.push_back(p);
        }

        collect(n, res);
        sort(res.begin(), res.end());
        return res;
    }

    /// Each prime factor of n is n itself or at least TRIAL_LIMIT, so below TRIAL_LIMIT^2 n is 1 or a prime
    void collect(long long n, vector<long long>& out){
        if (n == 1) return;

        if (n < TRIAL_LIMIT * TRIAL_LIMIT || prm::is_prime(n)){
            out.push_back(n);
            return;
        }

        long long d = n < (1LL << 61) ? split<true>(n) : split<false>(n);
        collect(d, out);
        collect(n / d, out);
    }

    /// A nontrivial factor of an odd composite n, LAZY needs n < 2^61
    template<bool LAZY>
    long long split(long long n){
        const int BATCH = 128;
        Montgomery mont(n);
        auto next = [&](unsigned long long x, unsigned long long c){
            if constexpr (LAZY){
                /// rho only needs x right modulo the unknown factor, so no conditional subtraction: 20% faster
                unsigned __int128 sq = (unsigned __int128)x * x;
                unsigned long long m = ((unsigned __int128)((unsigned long long)sq * mont.inv) * mont.n) >> 64;
                return (unsigned long long)(sq >> 64) + mont.n - m + c;  /// below 2^62 + n + c < 2^63 for x < 2^63
            }
            else return mont.add(mont.mul(x, x), c);
        };
        auto dist = [](unsigned long long x, unsigned long long y){ return x > y ? x - y : y - x; };

        while (true){
            unsigned long long c = rng() % (n - 1) + 1, y = rng() % n, x = y, ys = y, q = mont.to_mont(1), g = 1;

            for (long long r = 1; g == 1; r <<= 1){
                x = y;
                for (long long i = 0; i < r; i++) y = next(y, c);

                for (long long k = 0; k < r && g == 1; k += BATCH){
                    ys = y;
                    for (long long i = 0; i < BATCH && i < r - k; i++){
                        y = next(y, c);
                        q = mont.mul(q, dist(x, y));
                    }
                    g = gcd(q, (unsigned long long)n);  /// binary GCD in libstdc++: factorize 12-18% faster than __gcd
                }
            }

            if (g == (unsigned long long)n){  /// the batch overshot to a multiple of n, replay it one step at a time
                do{
                    ys = next(ys, c);
                    g = gcd(dist(x, ys), (unsigned long long)n);
                } while (g == 1);
            }

            if (g != (unsigned long long)n) return g;
        }
    }
};
/// END COPY pollard_rho

struct TwoSquares{
    using Gaussian = pair<long long, long long>;

    PollardRho rho;

    vector<pair<long long, long long>> pairs(long long n){
        if (n == 0) return {{0, 0}};

        /// One Gaussian integer of norm n per class of associates; unique factorization makes the classes distinct
        vector<Gaussian> zs = {{1, 0}};
        for (auto [p, e]: prime_powers(n)){
            if (p % 4 == 3 && e % 2 == 1) return {};

            vector<Gaussian> next;
            if (p == 2){
                for (auto z: zs){
                    for (int k = 0; k < e; k++) z = mul(z, {1, 1});
                    next.push_back(z);
                }
            }
            else if (p % 4 == 3){
                for (auto [a, b]: zs){
                    for (int k = 0; k < e / 2; k++) a *= p, b *= p;
                    next.push_back({a, b});
                }
            }
            else{
                /// p = pi conj(pi) with pi, conj(pi) not associates, so p^e has e + 1 classes pi^k conj(pi)^(e - k)
                Gaussian pi = prime_split(p);
                vector<Gaussian> pw = {{1, 0}};
                for (int k = 0; k < e; k++) pw.push_back(mul(pw.back(), pi));
                for (auto z: zs){
                    for (int k = 0; k <= e; k++) next.push_back(mul(z, mul(pw[k], {pw[e - k].first, -pw[e - k].second})));
                }
            }
            zs = next;
        }

        /// Rotating by i into x > 0, y >= 0 picks one of the 4 associates, so (0, x) comes from (x, 0)
        vector<pair<long long, long long>> res;
        for (auto [a, b]: zs){
            while (a <= 0 || b < 0) tie(a, b) = make_pair(-b, a);
            res.push_back({a, b});
            if (b == 0) res.push_back({0, a});
        }

        sort(res.begin(), res.end());
        return res;
    }

    long long r2(long long n){
        if (n == 0) return 1;

        long long res = 4;
        for (auto [p, e]: prime_powers(n)){
            if (p % 4 == 1) res *= e + 1;
            if (p % 4 == 3 && e % 2 == 1) return 0;
        }

        return res;
    }

    vector<pair<long long, int>> prime_powers(long long n){
        vector<pair<long long, int>> res;
        for (long long p: rho.factorize(n)){
            if (!res.empty() && res.back().first == p) res.back().second++;
            else res.push_back({p, 1});
        }
        return res;
    }

    /// Each product here has norm dividing n, so every component and partial product stays within sqrt(n)
    static Gaussian mul(Gaussian z, Gaussian w){
        return {z.first * w.first - z.second * w.second, z.first * w.second + z.second * w.first};
    }

    /// a + bi with a^2 + b^2 = p, for a prime p = 1 mod 4
    static Gaussian prime_split(long long p){
        Montgomery mont(p);
        unsigned long long minus_one = mont.n - mont.to_mont(1), x = 0;
        for (unsigned long long c = 2; mont.mul(x, x) != minus_one; c++) x = mont.pow(mont.to_mont(c), (p - 1) / 4);

        /// Brillhart 1972: for t^2 = -1 mod p, the first two Euclid remainders of (p, t) below sqrt(p) are a and b
        long long r0 = p, r1 = mont.reduce(x);
        while ((__int128)r1 * r1 > p){
            long long rem = r0 % r1;
            r0 = r1, r1 = rem;
        }
        return {r1, r0 % r1};
    }
};

int main(){
    /// 25 = 0 + 25 = 9 + 16, and 25 is a square so both (0, 5) and (5, 0) appear
    TwoSquares squares;
    vector<pair<long long, long long>> expected = {{0, 5}, {3, 4}, {4, 3}, {5, 0}};
    assert(squares.pairs(25) == expected);
    assert(squares.r2(25) == 12);  /// 4 points on the axes and 8 with x, y != 0

    /// 3 divides 21 once, so there are no pairs
    assert(squares.pairs(21).empty());
    assert(squares.r2(21) == 0);

    return 0;
}
