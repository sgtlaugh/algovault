#include "../common.h"

#define main library_main
#include "../../code_library/number_theory/chinese_remainder_theorem.cpp"
#undef main

const int64_t INF = numeric_limits<int64_t>::max();
const __int128 LIMIT = (__int128)1 << 63;

/// Pairwise coprime moduli whose product stays below limit
vector<int64_t> random_moduli(int count, int64_t max_mod, __int128 limit = LIMIT){
    vector<int64_t> mods;
    __int128 prod = 1;

    for (int tries = 0; (int)mods.size() < count && tries < 100; tries++){
        int64_t m = stress::rand_int(1, max_mod);
        bool coprime = true;
        for (auto x : mods) coprime &= __gcd(x, m) == 1;
        if (coprime && prod <= (limit - 1) / m) mods.push_back(m), prod *= m;
    }

    return mods;
}

/// Moduli sharing factors, lcm below 2^63: a common factor times small cofactors, or divisors of 2^63 - 1
vector<int64_t> random_shared_moduli(int count){
    const vector<int64_t> inf_factors = {7, 7, 73, 127, 337, 92737, 649657};

    while (true){
        vector<int64_t> mods;
        int64_t g = stress::rand_int(1, 1LL << stress::rand_int(0, 62));
        bool divisors = stress::rand_int(0, 3) == 0;

        for (int i = 0; i < count; i++){
            int64_t m = 1;
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
        for (auto m : mods) if (lcm < LIMIT) lcm = lcm / __gcd((int64_t)(lcm % m), m) * m;
        if (lcm < LIMIT) return mods;
    }
}

int64_t random_rem(){
    int kind = stress::rand_int(0, 3);
    if (kind == 0) return INT64_MIN;
    if (kind == 1) return INF;
    return stress::rand_int(INT64_MIN, INF);
}

bool solves(__int128 x, const vector<int64_t>& rems, const vector<int64_t>& mods){
    for (size_t i = 0; i < mods.size(); i++) if (((x - rems[i]) % mods[i] + mods[i]) % mods[i]) return false;
    return true;
}

/// Little-endian base 2^32 big integer, enough for the Garner reference
using Big = vector<uint32_t>;

int64_t big_mod(const Big& a, int64_t m){
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
        vector<int64_t> rems;
        for (auto m : mods) rems.push_back(stress::rand_int(-3 * m, 3 * m));

        int64_t expected = 0;
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
        vector<int64_t> rems;
        for (size_t i = 0; i < mods.size(); i++) rems.push_back(random_rem());

        int64_t x = CRT(rems, mods);
        assert(0 <= x && x < prod && solves(x, rems, mods));
    }
}

void tiny_crt_vs_scan(){
    /// Arbitrary moduli: scan [0, lcm) for the first solution, -1 when the scan finds none
    for (long long it = 0; it < stress::scaled(30000); it++){
        int count = stress::rand_int(1, 4);
        vector<int64_t> mods, rems;
        int64_t lcm = 1;
        while ((int)mods.size() < count){
            int64_t m = stress::rand_int(1, 16);
            if (lcm / __gcd(lcm, m) * m > 2000) break;
            mods.push_back(m), lcm = lcm / __gcd(lcm, m) * m;
        }

        int64_t planted = stress::rand_int(0, lcm - 1);
        for (auto m : mods) rems.push_back(it % 2 ? stress::rand_int(-3 * m, 3 * m) : planted + m * stress::rand_int(-3, 3));

        int64_t expected = 0;
        while (expected < lcm && !solves(expected, rems, mods)) expected++;
        if (expected == lcm) expected = -1;
        assert(CRT(rems, mods) == expected);
    }
}

void large_crt_vs_planted(){
    /// Shared factors up to the lcm limit: a planted solution must come back as planted % lcm, random residues follow the pairwise gcd rule
    for (long long it = 0; it < stress::scaled(100000); it++){
        auto mods = random_shared_moduli(stress::rand_int(1, 4));
        __int128 lcm = 1;
        for (auto m : mods) lcm = lcm / __gcd((int64_t)(lcm % m), m) * m;

        vector<int64_t> rems;
        if (it % 2 == 0){
            int64_t planted = stress::rand_int(INT64_MIN, INF);
            for (auto m : mods) rems.push_back(planted % m);
            int64_t expected = ((planted % lcm) + lcm) % lcm;
            assert(CRT(rems, mods) == expected);
            continue;
        }

        for (size_t i = 0; i < mods.size(); i++) rems.push_back(it % 4 == 1 ? random_rem() : stress::rand_int(-1, 1));
        bool consistent = true;
        for (size_t i = 0; i < mods.size(); i++){
            for (size_t j = 0; j < i; j++){
                int64_t g = __gcd(mods[i], mods[j]);
                consistent &= ((__int128)rems[i] - rems[j]) % g == 0;
            }
        }

        int64_t x = CRT(rems, mods);
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
        int64_t prod = 1;
        for (auto m : mods) prod *= m;
        vector<int64_t> rems;
        for (auto m : mods) rems.push_back(stress::rand_int(-3 * m, 3 * m));

        int64_t x = 0;
        while (!solves(x, rems, mods)) x++;
        int64_t mod = it % 3 ? stress::rand_int(1, 50) : stress::rand_int(1, INF);
        assert(garner(rems, mods, mod) == x % mod);
    }
}

void big_garner_vs_bignum(){
    /// Products far beyond 2^64: the reference builds x one congruence at a time by scanning the next digit
    for (long long it = 0; it < stress::scaled(5000); it++){
        auto mods = random_moduli(stress::rand_int(1, 25), 1000, (__int128)1 << 126);
        vector<int64_t> rems;
        for (size_t i = 0; i < mods.size(); i++) rems.push_back(random_rem());

        Big x = {0}, prod = {1};
        for (size_t i = 0; i < mods.size(); i++){
            int64_t m = mods[i], want = ((rems[i] % m) + m) % m;
            int64_t xm = big_mod(x, m), pm = big_mod(prod, m);
            uint32_t k = 0;
            while ((xm + k * pm) % m != want) k++;
            big_add_mul(x, prod, k);

            Big next = {0};
            big_add_mul(next, prod, m);
            prod = next;
        }

        int64_t mod = it % 3 == 0 ? INF : stress::rand_int(1, it % 3 == 1 ? 1000000 : INF);
        assert(garner(rems, mods, mod) == big_mod(x, mod));
    }
}

void huge_garner_properties(){
    /// Moduli near 2^63: x % mods[j] must give back rems[j], and the answer cannot depend on the order of the congruences
    for (long long it = 0; it < stress::scaled(20000); it++){
        int count = stress::rand_int(1, 5);
        vector<int64_t> mods, rems;
        while ((int)mods.size() < count){
            int64_t m = INF - stress::rand_int(0, 1000000);
            bool coprime = true;
            for (auto x : mods) coprime &= __gcd(x, m) == 1;
            if (coprime) mods.push_back(m);
        }

        for (size_t i = 0; i < mods.size(); i++) rems.push_back(random_rem());

        for (size_t j = 0; j < mods.size(); j++) assert(garner(rems, mods, mods[j]) == ((__int128)rems[j] % mods[j] + mods[j]) % mods[j]);

        int64_t mod = stress::rand_int(1, INF);
        auto rev_rems = rems, rev_mods = mods;
        reverse(rev_rems.begin(), rev_rems.end());
        reverse(rev_mods.begin(), rev_mods.end());
        assert(garner(rems, mods, mod) == garner(rev_rems, rev_mods, mod));
    }
}

int64_t pow_mod(int64_t b, int64_t e, int64_t m){
    int64_t r = 1 % m;
    for (b %= m; e; e >>= 1, b = (__int128)b * b % m) if (e & 1) r = (__int128)r * b % m;
    return r;
}

void ntt_garner_recombination(){
    /// One 2^21 three-prime NTT convolution, one Garner reused for every coefficient against the textbook CRT sum
    const vector<int64_t> primes = {998244353, 167772161, 469762049};
    const int64_t mod = 1000000007;
    const int n = 1 << 21;

    __int128 prod = (__int128)primes[0] * primes[1] * primes[2];
    vector<__int128> basis;
    for (auto p : primes){
        __int128 rest = prod / p;
        basis.push_back(rest * pow_mod(rest % p, p - 2, p) % prod);
    }

    Garner recombine(primes, mod);
    vector<int64_t> rems(3);
    for (int it = 0; it < n; it++){
        __int128 x = 0;
        for (int i = 0; i < 3; i++){
            rems[i] = stress::rand_int(0, primes[i] - 1);
            x = (x + basis[i] * rems[i]) % prod;
        }
        assert(recombine(rems) == x % mod);
    }
}

int main(){
    tiny_coprime_vs_scan();
    large_coprime_congruences();
    tiny_crt_vs_scan();
    large_crt_vs_planted();
    tiny_garner_vs_scan();
    big_garner_vs_bignum();
    huge_garner_properties();
    ntt_garner_recombination();

    return 0;
}
