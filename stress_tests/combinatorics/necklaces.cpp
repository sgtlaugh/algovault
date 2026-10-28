#include "../common.h"

#define main bignum_main
#include "../../code_library/algebra/bignum.cpp"
#undef main

#define main library_main
#include "../../code_library/combinatorics/necklaces.cpp"
#undef main

typedef unsigned long long ull;

/// Orbits of all k^n colorings under rotation, and under rotation and reflection, each counted by its smallest member
pair<long long, long long> brute(int n, int k){
    long long total = 1;
    for (int i = 0; i < n; i++) total *= k;

    long long neck = 0, brace = 0;
    vector<int> c(n), t(n);
    for (long long code = 0; code < total; code++){
        long long x = code;
        for (int i = 0; i < n; i++) c[i] = x % k, x /= k;

        bool min_rot = true, min_all = true;
        for (int s = 0; s < n; s++){
            for (int dir = 0; dir < 2; dir++){
                for (int i = 0; i < n; i++) t[i] = dir ? c[(s - i + n) % n] : c[(s + i) % n];
                if (t < c){
                    min_all = false;
                    if (!dir) min_rot = false;
                }
            }
        }

        neck += min_rot, brace += min_all;
    }

    return {neck, brace};
}

ull mul_mod(ull a, ull b, ull m){
    return (unsigned __int128)a * b % m;
}

ull pow_mod(ull b, ull e, ull m){
    ull res = 1 % m;
    for (b %= m; e; e >>= 1, b = mul_mod(b, b, m)){
        if (e & 1) res = mul_mod(res, b, m);
    }
    return res;
}

/// Necklaces and bracelets modulo a prime m coprime to 2 n, from n's divisors found by plain trial division
pair<ull, ull> reference_mod(long long n, long long k, ull m){
    vector<long long> primes;
    long long rest = n;
    for (long long p = 2; p * p <= rest; p++){
        if (rest % p == 0) primes.push_back(p);
        while (rest % p == 0) rest /= p;
    }
    if (rest > 1) primes.push_back(rest);

    ull sum = 0;
    for (long long d = 1; d * d <= n; d++){
        if (n % d) continue;
        for (long long e : {d, n / d}){
            long long phi = e;
            for (long long p : primes) if (e % p == 0) phi = phi / p * (p - 1);
            sum = (sum + mul_mod(phi % m, pow_mod(k, n / e, m), m)) % m;
            if (d * d == n) break;
        }
    }

    ull refl = n % 2 ? mul_mod(n % m, pow_mod(k, (n + 1) / 2, m), m)
                     : mul_mod((n / 2) % m, (pow_mod(k, n / 2, m) + pow_mod(k, n / 2 + 1, m)) % m, m);
    ull neck = mul_mod(sum, pow_mod(n, m - 2, m), m);
    ull brace = mul_mod((sum + refl) % m, pow_mod(2 * n % m, m - 2, m), m);
    return {neck, brace};
}

void check_mod(long long n, long long k, ull m, ull neck, ull brace){
    if ((unsigned __int128)n * m <= ULLONG_MAX) assert((ull)necklaces(n, k, m) == neck);
    if ((unsigned __int128)2 * n * m <= ULLONG_MAX) assert((ull)bracelets(n, k, m) == brace);
}

bool is_prime(ull x){
    if (x < 2) return false;
    for (ull d = 2; d * d <= x; d++){
        if (x % d == 0) return false;
    }
    return true;
}

int main(){
    stress::expect_assert("combinatorics/necklaces.cpp", []{ necklaces<long long>(0, 2); });
    stress::expect_assert("combinatorics/necklaces.cpp", []{ bracelets<long long>(3, -1); });
    stress::expect_assert("combinatorics/necklaces.cpp", []{ necklaces(1LL << 32, 2, 1LL << 32); });
    stress::expect_assert("combinatorics/necklaces.cpp", []{ bracelets(1LL << 31, 2, 1LL << 32); });

    /// Every small ring against the brute force in every result type, and modulo every m <= 40, divisors of n included
    for (int n = 1; n <= 12; n++){
        for (int k = 0; k <= 5; k++){
            if (pow(k, n) * 2 * n * n > 1e7) continue;
            auto [neck, brace] = brute(n, k);
            assert(necklaces<long long>(n, k) == neck);
            assert(bracelets<long long>(n, k) == brace);
            assert(necklaces<__int128>(n, k) == neck);
            assert(bracelets<__int128>(n, k) == brace);
            assert(necklaces<Bignum>(n, k) == neck);
            assert(bracelets<Bignum>(n, k) == brace);
            for (long long m = 1; m <= 40; m++){
                assert(necklaces(n, k, m) == neck % m);
                assert(bracelets(n, k, m) == brace % m);
            }
        }
    }

    /// __int128 and Bignum agree modulo two coprime m near 10^18 while k^n < 2^115, so the counts are pinned exactly
    for (int n = 1; n <= 115; n++){
        for (int k : {2, 3, 7, 10}){
            if (n * log2(k) > 115) continue;
            for (long long m : {1000000000000000000LL, 999999999999999989LL}){
                assert(necklaces<Bignum>(n, k) % m == (long long)(necklaces<__int128>(n, k) % m));
                assert(bracelets<Bignum>(n, k) % m == (long long)(bracelets<__int128>(n, k) % m));
            }
        }
    }

    /// n times the answer right below 2^63, where squaring past k^n in the power would overflow long long
    assert(necklaces<long long>(2, 3037000499LL) == necklaces<Bignum>(2, 3037000499LL));
    assert(necklaces<long long>(62, 2) == necklaces<Bignum>(62, 2));
    assert(bracelets<long long>(62, 2) == bracelets<Bignum>(62, 2));

    /// Exact counts far past 128 bits reduced by m match the modular calls, which never see the big numbers
    for (int it = 0; it < stress::scaled(30); it++){
        long long n = stress::rand_int(1, 300), k = stress::rand_int(0, 1000), m = stress::rand_int(1, 1000000);
        assert(necklaces<Bignum>(n, k) % m == necklaces(n, k, m));
        assert(bracelets<Bignum>(n, k) % m == bracelets(n, k, m));
    }

    /// UVa 10294 bounds, n <= 50 and t <= 10, in long long wherever the answer has at most 11 digits
    for (int n = 1; n <= 50; n++){
        for (int t = 1; t <= 10; t++){
            Bignum neck = necklaces<Bignum>(n, t), brace = bracelets<Bignum>(n, t);
            if (neck >= 100000000000LL) continue;
            assert(necklaces<long long>(n, t) == neck);
            assert(bracelets<long long>(n, t) == brace);
        }
    }

    /// Random n against the independent trial division reference, for primes m coprime to 2 n
    vector<ull> big_primes = {1000000007, 998244353, 4294967291ULL, 4294967279ULL};
    for (int it = 0; it < stress::scaled(200); it++){
        long long n = it % 10 ? stress::rand_int(1, 1000000) : stress::rand_int(1, 1000000000000LL);
        long long k = stress::rand_int(0, 1000000000000000000LL);
        ull m = it % 3 ? big_primes[stress::rand_int(0, big_primes.size() - 1)] : 0;
        while (m == 0 || !is_prime(m)) m = stress::rand_int(3, 20000000);
        if (n % m == 0) continue;
        auto [neck, brace] = reference_mod(n, k, m);
        check_mod(n, k, m, neck, brace);
    }

    /// Highly composite n, the most divisors below 10^12, and products n m right below 2^64
    for (long long n : {963761198400LL, 735134400LL, 4294967291LL, 4294967296LL, 999999999989LL}){
        for (ull m : {3ULL, 1000003ULL, 9570181ULL, 19140293ULL, 4294967279ULL, 4294967291ULL}){
            if (n % m == 0 || (unsigned __int128)n * m > ULLONG_MAX) continue;
            auto [neck, brace] = reference_mod(n, 123456789, m);
            check_mod(n, 123456789, m, neck, brace);
        }
    }

    /// Prime n: every non-constant coloring has n distinct rotations, so N = (k^n - k) / n + k
    for (long long n : {2LL, 3LL, 13LL, 61LL, 1000003LL, 999999999989LL}){
        for (long long m : {2LL, 6LL, 1000003LL}){
            long long k = 5;
            ull mod = n * m, kn = pow_mod(k, n, mod);
            assert(necklaces(n, k, m) == (long long)((kn + mod - k % mod) % mod / n + k) % m);
        }
    }

    return 0;
}
