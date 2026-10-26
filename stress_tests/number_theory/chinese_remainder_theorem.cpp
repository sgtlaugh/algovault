#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/chinese_remainder_theorem.cpp"
#undef main

const long long INF = numeric_limits<long long>::max();
const __int128 LIMIT = (__int128)1 << 63;

/// Pairwise coprime moduli whose product stays below limit
vector<long long> random_moduli(int count, long long max_mod, __int128 limit = LIMIT){
    vector<long long> mods;
    __int128 prod = 1;

    for (int tries = 0; (int)mods.size() < count && tries < 100; tries++){
        long long m = stress::rand_int(1, max_mod);
        bool coprime = true;
        for (auto x : mods) coprime &= __gcd(x, m) == 1;
        if (coprime && prod <= (limit - 1) / m) mods.push_back(m), prod *= m;
    }

    return mods;
}

/// Moduli sharing factors, lcm below 2^63: a common factor times small cofactors, or divisors of 2^63 - 1
vector<long long> random_shared_moduli(int count){
    const vector<long long> inf_factors = {7, 7, 73, 127, 337, 92737, 649657};

    while (true){
        vector<long long> mods;
        long long g = stress::rand_int(1, 1LL << stress::rand_int(0, 62));
        bool divisors = stress::rand_int(0, 3) == 0;

        for (int i = 0; i < count; i++){
            long long m = 1;
            if (divisors){
                for (auto p : inf_factors) if (stress::rand_int(0, 1)) m *= p;
            }
            else{
                __int128 v = (__int128)g * stress::rand_int(1, 1LL << stress::rand_int(0, 20));
                if (v >= LIMIT) break;
                m = v;
            }
            mods.push_back(m);
        }

        if ((int)mods.size() < count) continue;

        __int128 lcm = 1;
        for (auto m : mods) if (lcm < LIMIT) lcm = lcm / __gcd((long long)(lcm % m), m) * m;
        if (lcm < LIMIT) return mods;
    }
}

long long random_rem(){
    int kind = stress::rand_int(0, 3);
    if (kind == 0) return INT64_MIN;
    if (kind == 1) return INF;
    return stress::rand_int(INT64_MIN, INF);
}

bool solves(__int128 x, const vector<long long>& rems, const vector<long long>& mods){
    for (size_t i = 0; i < mods.size(); i++) if (((x - rems[i]) % mods[i] + mods[i]) % mods[i]) return false;
    return true;
}

/// Little-endian base 2^32 big integer, enough for the Garner reference
using Big = vector<uint32_t>;

long long big_mod(const Big& a, long long m){
    unsigned __int128 rem = 0;
    for (int i = (int)a.size() - 1; i >= 0; i--) rem = ((rem << 32) | a[i]) % m;
    return rem;
}

/// a += b * k
void big_add_mul(Big& a, const Big& b, uint32_t k){
    uint64_t carry = 0;
    if (a.size() < b.size()) a.resize(b.size(), 0);

    for (size_t i = 0; i < a.size(); i++){
        uint64_t cur = a[i] + carry + (i < b.size() ? (uint64_t)b[i] * k : 0);
        a[i] = cur, carry = cur >> 32;
    }

    if (carry) a.push_back(carry);
}

void tiny_coprime_vs_scan(){
    /// Tiny co-prime moduli: the unique answer in [0, prod) is found by brute force
    for (long long it = 0; it < stress::scaled(20000); it++){
        auto mods = random_moduli(stress::rand_int(1, 4), 12);
        vector<long long> rems;
        for (auto m : mods) rems.push_back(stress::rand_int(-3 * m, 3 * m));

        long long expected = 0;
        while (!solves(expected, rems, mods)) expected++;
        assert(CRT(rems, mods) == expected);
    }
}

void large_coprime_congruences(){
    /// Large co-prime moduli: check the congruences exactly
    for (long long it = 0; it < stress::scaled(100000); it++){
        auto mods = random_moduli(stress::rand_int(1, 3), it % 4 == 0 ? INF : it % 2 ? 3000000000LL : 1000000);
        __int128 prod = 1;
        for (auto m : mods) prod *= m;
        vector<long long> rems;
        for (size_t i = 0; i < mods.size(); i++) rems.push_back(random_rem());

        long long x = CRT(rems, mods);
        assert(0 <= x && x < prod && solves(x, rems, mods));
    }
}

void tiny_crt_vs_scan(){
    /// Arbitrary moduli: scan [0, lcm) for the first solution, -1 when the scan finds none
    for (long long it = 0; it < stress::scaled(30000); it++){
        int count = stress::rand_int(1, 4);
        vector<long long> mods, rems;
        long long lcm = 1;
        while ((int)mods.size() < count){
            long long m = stress::rand_int(1, 16);
            if (lcm / __gcd(lcm, m) * m > 2000) break;
            mods.push_back(m), lcm = lcm / __gcd(lcm, m) * m;
        }

        long long planted = stress::rand_int(0, lcm - 1);
        for (auto m : mods) rems.push_back(it % 2 ? stress::rand_int(-3 * m, 3 * m) : planted + m * stress::rand_int(-3, 3));

        long long expected = 0;
        while (expected < lcm && !solves(expected, rems, mods)) expected++;
        if (expected == lcm) expected = -1;
        assert(CRT(rems, mods) == expected);
    }
}

void large_crt_vs_planted(){
    /// Shared factors up to the lcm limit: a planted x comes back as x % lcm, random residues follow the gcd rule
    for (long long it = 0; it < stress::scaled(100000); it++){
        auto mods = random_shared_moduli(stress::rand_int(1, 4));
        __int128 lcm = 1;
        for (auto m : mods) lcm = lcm / __gcd((long long)(lcm % m), m) * m;

        vector<long long> rems;
        if (it % 2 == 0){
            long long planted = stress::rand_int(INT64_MIN, INF);
            for (auto m : mods) rems.push_back(planted % m);
            long long expected = ((planted % lcm) + lcm) % lcm;
            assert(CRT(rems, mods) == expected);
            continue;
        }

        for (size_t i = 0; i < mods.size(); i++) rems.push_back(it % 4 == 1 ? random_rem() : stress::rand_int(-1, 1));
        bool consistent = true;
        for (size_t i = 0; i < mods.size(); i++){
            for (size_t j = 0; j < i; j++){
                long long g = __gcd(mods[i], mods[j]);
                consistent &= ((__int128)rems[i] - rems[j]) % g == 0;
            }
        }

        long long x = CRT(rems, mods);
        if (!consistent){
            assert(x == -1);
            continue;
        }
        assert(0 <= x && x < lcm && solves(x, rems, mods));
    }
}

void tiny_garner_vs_scan(){
    for (long long it = 0; it < stress::scaled(10000); it++){
        auto mods = random_moduli(stress::rand_int(0, 3), 30);
        long long prod = 1;
        for (auto m : mods) prod *= m;
        vector<long long> rems;
        for (auto m : mods) rems.push_back(stress::rand_int(-3 * m, 3 * m));

        long long x = 0;
        while (!solves(x, rems, mods)) x++;
        long long mod = it % 3 ? stress::rand_int(1, 50) : stress::rand_int(1, INF);
        assert(Garner(mods, mod)(rems) == x % mod);
    }
}

void big_garner_vs_bignum(){
    /// Products far beyond 2^64: the reference builds x one congruence at a time by scanning the next digit
    for (long long it = 0; it < stress::scaled(5000); it++){
        auto mods = random_moduli(stress::rand_int(1, 25), 1000, (__int128)1 << 126);
        vector<long long> rems;
        for (size_t i = 0; i < mods.size(); i++) rems.push_back(random_rem());

        Big x = {0}, prod = {1};
        for (size_t i = 0; i < mods.size(); i++){
            long long m = mods[i], want = ((rems[i] % m) + m) % m;
            long long xm = big_mod(x, m), pm = big_mod(prod, m);
            uint32_t k = 0;
            while ((xm + k * pm) % m != want) k++;
            big_add_mul(x, prod, k);

            Big next = {0};
            big_add_mul(next, prod, m);
            prod = next;
        }

        long long mod = it % 3 == 0 ? INF : stress::rand_int(1, it % 3 == 1 ? 1000000 : INF);
        assert(Garner(mods, mod)(rems) == big_mod(x, mod));
    }
}

void huge_garner_properties(){
    /// Moduli near 2^63: x % mods[j] must give back rems[j], and the answer must not depend on the congruence order
    for (long long it = 0; it < stress::scaled(20000); it++){
        int count = stress::rand_int(1, 5);
        vector<long long> mods, rems;
        while ((int)mods.size() < count){
            long long m = INF - stress::rand_int(0, 1000000);
            bool coprime = true;
            for (auto x : mods) coprime &= __gcd(x, m) == 1;
            if (coprime) mods.push_back(m);
        }

        for (size_t i = 0; i < mods.size(); i++) rems.push_back(random_rem());

        for (size_t j = 0; j < mods.size(); j++) assert(Garner(mods, mods[j])(rems) == ((__int128)rems[j] % mods[j] + mods[j]) % mods[j]);

        long long mod = stress::rand_int(1, INF);
        auto rev_rems = rems, rev_mods = mods;
        reverse(rev_rems.begin(), rev_rems.end());
        reverse(rev_mods.begin(), rev_mods.end());
        assert(Garner(mods, mod)(rems) == Garner(rev_mods, mod)(rev_rems));
    }
}

long long pow_mod(long long b, long long e, long long m){
    long long r = 1 % m;
    for (b %= m; e; e >>= 1, b = (__int128)b * b % m) if (e & 1) r = (__int128)r * b % m;
    return r;
}

void ntt_garner_recombination(){
    /// One 2^21 three-prime NTT convolution, one Garner reused for every coefficient against the textbook CRT sum
    const vector<long long> primes = {998244353, 167772161, 469762049};
    const long long mod = 1000000007;
    const int n = 1 << 21;

    __int128 prod = (__int128)primes[0] * primes[1] * primes[2];
    vector<__int128> basis;
    for (auto p : primes){
        __int128 rest = prod / p;
        basis.push_back(rest * pow_mod(rest % p, p - 2, p) % prod);
    }

    Garner recombine(primes, mod);
    vector<long long> rems(3);
    for (int it = 0; it < n; it++){
        __int128 x = 0;
        for (int i = 0; i < 3; i++){
            rems[i] = stress::rand_int(0, primes[i] - 1);
            x = (x + basis[i] * rems[i]) % prod;
        }
        assert(recombine(rems) == x % mod);
    }
}

void edge_cases(){
    /// Empty systems, mod 1, residues at INT64_MIN and INF, and lcm or target mods right below 2^63
    assert(CRT({}, {}) == 0);
    assert(CRT({5}, {1}) == 0);
    assert(CRT({-1, -1}, {4, 6}) == 11);
    assert(CRT({7, 2}, {5, 5}) == 2);
    assert(CRT({1, 2}, {5, 5}) == -1);
    assert(CRT({1, 2, 3}, {2, 3, 4}) == 11);
    assert(CRT({0, 2, 3}, {2, 3, 4}) == -1);
    assert(CRT({5, 3 * (1LL << 60) + 5}, {3 * (1LL << 61), 3 * (1LL << 60)}) == 5);
    assert(CRT({5, 6}, {3 * (1LL << 61), 3 * (1LL << 60)}) == -1);
    assert(CRT({-51146671749066063, 9680838726321}, {1099511627776LL * 59049, 1073741824LL * 1701}) == 78703452468023985);
    assert(CRT({-1}, {INF}) == INF - 1);
    assert(CRT({INT64_MIN}, {INF}) == INF - 1);
    assert(CRT({0, 1}, {153092023, 60247241209}) == 2695413193155969174);
    assert(CRT({INT64_MIN, INF}, {153092023, 60247241209}) == 2695413193155969173);

    assert(Garner({}, 10)({}) == 0);
    assert(Garner({3, 5, 7}, 1)({2, 3, 2}) == 0);
    assert(Garner({1, 1, 7}, 10)({5, 6, 3}) == 3);
    assert(Garner({1000000007, 998244353, 1000000009}, 1000000007)({1, 2, 3}) == 1);
    assert(Garner({1000000007, 998244353, 1000000009}, (1LL << 61) - 1)({1, 2, 3}) == 1131585747248034387);
    assert(Garner({1000000007, 998244353, 1000000009}, INF)({1, 2, 3}) == 1131585747206765391);

    vector<long long> big_rems = {-1, 100000000000000000, 1LL << 62, INT64_MIN};
    vector<long long> big_mods = {1000000000000000003, 1000000000000000009, (1LL << 61) - 1, INF};
    assert(Garner(big_mods, INF)(big_rems) == INF - 1);
    assert(Garner(big_mods, 1000000007)(big_rems) == 746669198);
    assert(Garner(big_mods, 2)(big_rems) == 0);
}

void rejected_inputs(){
    /// Garner on moduli sharing a factor and CRT past an lcm of 2^63 must assert, not return a wrong answer
    const string lib = "number_theory/chinese_remainder_theorem.cpp";
    stress::expect_assert(lib, []{ Garner({4, 6}, 10); });
    stress::expect_assert(lib, []{ Garner({7, 3, 21}, 10); });
    stress::expect_assert(lib, []{ CRT({1, 2}, {4000000007LL, 4000000009LL}); });
    stress::expect_assert(lib, []{ CRT({0, 0, 0}, {1LL << 62, 3, 2}); });
}

void garner_wide_sums(){
    /// Digits and coefficients near 2^63 for 12 moduli, so a prefix sum reduced too lazily overflows 128 bits
    const long long target = 9223372036854775783;
    for (long long it = 0; it < stress::scaled(100); it++){
        vector<long long> mods = {target - stress::rand_int(1, 1LL << 40)};
        long long coef = mods[0], coef_inv = pow_mod(coef, target - 2, target);
        while (mods.size() < 12){
            long long next = target - stress::rand_int(1, 1LL << 40);
            long long m = (__int128)next * coef_inv % target;
            bool usable = m > target - target / 16;
            for (auto x : mods) usable = usable && __gcd(x, m) == 1;
            if (usable) mods.push_back(m), coef = next, coef_inv = pow_mod(coef, target - 2, target);
        }

        vector<long long> rems(mods.size(), -1);  /// x = product - 1, every digit at its maximum m - 1
        assert(Garner(mods, target)(rems) == coef - 1);
    }
}

int main(){
    edge_cases();
    rejected_inputs();
    tiny_coprime_vs_scan();
    large_coprime_congruences();
    tiny_crt_vs_scan();
    large_crt_vs_planted();
    tiny_garner_vs_scan();
    big_garner_vs_bignum();
    huge_garner_properties();
    garner_wide_sums();
    ntt_garner_recombination();

    return 0;
}
