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
    vector<uint32_t> dp;

    uint64_t get_tz(uint64_t x){
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
            if (d & 1) res = res * dp[mod - 1] % mod;
            res = res * dp[x - d * mod] % mod;
            x /= p;
        }

        return res * dp[x] % mod;
    }

    BinomialPrimePow(uint64_t p, uint32_t q, int64_t max_n): p(p), q(q){
        /// p^q divides the 64-bit mod, so this never overflows before the assert
        uint64_t pw = 1;
        for (uint32_t i = 0; i < q; i++) pw *= p;
        assert(pw < (1ULL << 32));
        this->mod = pw;
        /// While n <= max_n < p^q, fact() only ever reads dp[0..max_n], so a longer table would never be touched
        this->len = min<uint64_t>(pw, (uint64_t)max_n + 1);
    }

    int64_t binomial(int64_t n, int64_t k){
        uint64_t t = get_tz(n) - get_tz(k) - get_tz(n - k);
        if (t >= q) return 0;

        if (!dp.size()){
            dp.resize(len, 1);
            for (uint32_t i = 1, r = 1; i < len; i++, r++){
                if (r == p) r = 0;
                dp[i] = (uint64_t)dp[i - 1] * (r ? i : 1) % mod;
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
    vector<BinomialPrimePow> dp;
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
            uint64_t x = (uint64_t)rem[i] * invs[i] % dp[i].mod * mods[i];
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
            if (c) dp.push_back(BinomialPrimePow(i, c, max_n));
        }
        if (m > 1) dp.push_back(BinomialPrimePow(m, 1, max_n));

        for (auto d: dp){
            mods.push_back(mod / d.mod);
            invs.push_back(mod_inverse(mod / d.mod % d.mod, d.mod));
        }
    }

    int64_t binomial(int64_t n, int64_t k){
        if (k > n || n < 0 || k < 0 || !dp.size()) return 0;
        assert(n <= max_n);
        if (n == k || k == 0) return 1;

        vector<uint32_t> rem;
        for (uint32_t i = 0; i < dp.size(); i++){
            rem.push_back(dp[i].binomial(n, k));
        }

        return chinese_remainder(rem);
    }
};

int main(){
    Binomial bin(1000000007, 31);
    long long C[32][32];

    for (int i = 0; i < 32; i++){
        for (int k = 0; k < 32; k++) C[i][k] = bin.binomial(i, k);
    }

    long long n;
    while (scanf("%lld", &n) == 1){
        int k = 0;
        while (n >= C[31][k]) n -= C[31][k++];

        long long res = 0;
        for (int i = 30; i >= 0; i--){
            res <<= 1;
            if (n >= C[i][k]){
                res |= 1;
                n -= C[i][k--];
            }
        }
        printf("%lld\n", res);
    }

    return 0;
}
