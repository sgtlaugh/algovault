/***
 *
 * Binomial Coefficients
 * C(n, k) modulo any mod, prime or not, by Lucas per prime power of mod joined with the Chinese remainder theorem
 *
 * Complexity, for mod = p1^q1 * ... * pm^qm and table sizes Pi = min(pi^qi, max_n + 1):
 *   - constructor: O(sqrt(mod)) trial division, instant below 2^32, about 50 s for two prime factors near 2^32
 *   - tables: O(Pi) time and 4 Pi bytes per prime power, built by the first query that needs them
 *   - query: O(m (log n + log mod))
 *
 * Binomial bin(mod, max_n = INT64_MAX): 1 <= mod < 2^64, every prime power below 2^32, queries need n <= max_n
 * bin.binomial(n, k): C(n, k) mod mod for n < 2^63, 0 when n < 0, k < 0 or k > n, always 0 for mod 1
 *
 * Without max_n a prime mod near 2e9 builds an 8 GB table, so pass max_n whenever n stays small
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct BinomialPrimePow{
    uint32_t p, q, mod, len;
    /// free_fact[i] = the product of 1..i skipping multiples of p, modulo p^q
    vector<uint32_t> free_fact;

    uint64_t legendre(uint64_t x){
        uint64_t res = 0;
        while (x){
            x /= p;
            res += x;
        }

        return res;
    }

    uint64_t expo(uint64_t x, uint64_t n){
        x %= mod;
        uint64_t res = 1;

        while (n){
            if (n & 1) res = res * x % mod;
            x = x * x % mod;
            n >>= 1;
        }

        return res;
    }

    uint64_t fact(uint64_t x){
        uint64_t d, res = 1;

        while (x >= p){
            d = x / mod;
            if (d & 1) res = res * free_fact[mod - 1] % mod;
            res = res * free_fact[x - d * mod] % mod;
            x /= p;
        }

        return res * free_fact[x] % mod;
    }

    BinomialPrimePow(uint64_t p, uint32_t q, int64_t max_n): p(p), q(q){
        /// p^q divides the 64-bit mod, so this never overflows before the assert
        uint64_t pw = 1;
        for (uint32_t i = 0; i < q; i++) pw *= p;
        assert(pw < (1ULL << 32));
        this->mod = pw;
        /// While n <= max_n < p^q, fact() only ever reads free_fact[0..max_n], so a longer table would never be touched
        this->len = min<uint64_t>(pw, (uint64_t)max_n + 1);
    }

    int64_t binomial(int64_t n, int64_t k){
        uint64_t t = legendre(n) - legendre(k) - legendre(n - k);
        if (t >= q) return 0;

        if (!free_fact.size()){
            free_fact.resize(len, 1);
            for (uint32_t i = 1, r = 1; i < len; i++, r++){
                if (r == p) r = 0;
                free_fact[i] = (uint64_t)free_fact[i - 1] * (r ? i : 1) % mod;
            }
        }

        uint64_t m = fact(k) * fact(n - k) % mod;
        uint64_t v = fact(n) * expo(m, (mod / p) * (p - 1) - 1) % mod;
        return v * expo(p, t) % mod;
    }
};

struct Binomial{
    uint64_t mod;
    int64_t max_n;
    vector<BinomialPrimePow> prime_powers;
    vector<uint64_t> mods;
    vector<uint32_t> invs;

    uint32_t extended_gcd(uint32_t a, uint32_t b, uint32_t& x, uint32_t& y){
        if (!b){
            y = 0, x = 1;
            return a;
        }

        uint32_t g = extended_gcd(b, a % b, y, x);
        y -= ((a / b) * x);
        return g;
    }

    uint32_t mod_inverse(uint32_t a, uint32_t m){
        uint32_t x, y, inv;
        extended_gcd(a, m, x, y);
        inv = (int32_t)x < 0 ? x + m : x;  /// |x| <= m / 2 < 2^31, (x + m) % m wraps in 32 bits for m near 2^32
        return inv;
    }

    uint64_t chinese_remainder(const vector<uint32_t>& rem){
        uint64_t res = 0;
        for (uint32_t i = 0; i < rem.size(); i++){
            uint64_t x = (uint64_t)rem[i] * invs[i] % prime_powers[i].mod * mods[i];
            res = res >= mod - x ? res - (mod - x) : res + x;
        }

        return res;
    }

    Binomial(uint64_t mod, int64_t max_n = INT64_MAX): mod(mod), max_n(max_n){
        assert(max_n >= 0);
        uint64_t m = mod;

        for (uint64_t i = 2; i <= m / i; i++){
            uint32_t c = 0;
            while (m % i == 0){
                c++;
                m /= i;
            }
            if (c) prime_powers.push_back(BinomialPrimePow(i, c, max_n));
        }
        if (m > 1) prime_powers.push_back(BinomialPrimePow(m, 1, max_n));

        for (const auto& pp: prime_powers){
            mods.push_back(mod / pp.mod);
            invs.push_back(mod_inverse(mod / pp.mod % pp.mod, pp.mod));
        }
    }

    int64_t binomial(int64_t n, int64_t k){
        if (k > n || n < 0 || k < 0 || !prime_powers.size()) return 0;
        assert(n <= max_n);
        if (n == k || k == 0) return 1;

        vector<uint32_t> rem;
        for (uint32_t i = 0; i < prime_powers.size(); i++){
            rem.push_back(prime_powers[i].binomial(n, k));
        }
        return chinese_remainder(rem);
    }
};

int main(){
    /// Every C(n, k) with n, k < 30 against Pascal's triangle, for each modulus below 30, prime or not
    int64_t pascal[30][30];
    for (int mod = 1; mod < 30; mod++){
        memset(pascal, 0, sizeof(pascal));
        Binomial bin(mod);

        for (int n = 0; n < 30; n++){
            for (int k = 0; k <= n; k++){
                if (k == 0 || n == k) pascal[n][k] = 1 % mod;
                else pascal[n][k] = (pascal[n - 1][k - 1] + pascal[n - 1][k]) % mod;
            }
        }

        for (int n = 0; n < 30; n++){
            for (int k = 0; k < 30; k++){
                assert(pascal[n][k] == bin.binomial(n, k));
            }
        }
    }

    assert(Binomial(253).binomial(10, 5) == 252);
    assert(Binomial(997 * 1000003).binomial(123456789012345LL, 1000000009) == 140000420);
    assert(Binomial(1000003).binomial(1000000000000000000LL, 1000000000) == 0);
    assert(Binomial(666666667).binomial(1000000000000000000LL, 1000000000) == 0);
    assert(Binomial(1000000007 * 997LL).binomial(100000000000000000LL, 10000000009LL) == 353000002471LL);
    assert(Binomial(1000000007 * 997LL).binomial(1000000000000000000LL, 10000000009LL) == 361000002527LL);
    assert(Binomial(1000000007 * 997LL).binomial(1000000000000000000LL, 1000) == 625000004375LL);

    /// 4294967291 is the largest prime below 2^32, C(p, 1) = p and p < 3p
    assert(Binomial(3 * 4294967291ULL).binomial(4294967291LL, 1) == 4294967291LL);

    Binomial wide(3 * 4294967291ULL);
    for (size_t i = 0; i < wide.prime_powers.size(); i++) assert(wide.mods[i] % wide.prime_powers[i].mod * wide.invs[i] % wide.prime_powers[i].mod == 1);

    /// 2^63 - 1 is all ones in binary, so by Lucas every C(2^63 - 1, k) is odd
    assert(Binomial(2).binomial(INT64_MAX, 123456789) == 1);

    /// Uncapped, these would build tables of 2e9, 1.16e9 and 1e9 entries, max_n keeps each at 100001
    assert(Binomial(1999999973, 100000).binomial(100000, 50000) == 1615655902);
    assert(Binomial(1162261467, 100000).binomial(100000, 33333) == 641856042);  /// 3^19
    assert(Binomial(2 * 999999937, 100000).binomial(100000, 33333) == 166212410);

    return 0;
}
