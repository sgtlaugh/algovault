/***
 *
 * Binomial co-efficients (n choose k) modulo mod
 *
 * Generalization of the Lucas theorem, which calculates binomial mod p^q
 *
 * Let mod = p1^q1 * p2^q2 * ... * pm^qm
 *
 * Then, it will take O(sqrt(mod)) + O(max(p_i^q_i)) to pre-process
 * And answer each n choose k query in O(sum of log(k) * log(p_i^q_i))
 * Or roughly O(m * log(k) * log(mod))
 *
***/

#include <bits/stdc++.h>

using namespace std;

struct BinomialPP{
    uint32_t p, q, mod;
    vector<uint32_t> dp;

    BinomialPP(uint32_t p, uint32_t q): p(p), q(q){
        this->mod = pow(p, q) + 0.25;
    }

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

        return res % mod;
    }

    uint64_t fact(uint64_t x){
        uint64_t res = expo(dp[mod - 1], x / mod);
        if (x >= p) res = res * fact(x / p) % mod;
        return res * dp[x % mod] % mod;
    }

    int64_t binomial(int64_t n, int64_t k){
        uint64_t t = get_tz(n) - get_tz(k) - get_tz(n - k);
        if (t >= q) return 0;

        if (!dp.size()){
            dp.resize(mod, 1);
            for (uint32_t i = 1, r = 1; i < mod; i++, r++){
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
    vector<BinomialPP> dp;

    Binomial(uint64_t mod){
        for (uint32_t i = 2; i * i <= mod; i++){
            uint32_t c = 0;
            while (mod % i == 0){
                c++;
                mod /= i;
            }
            if (c) dp.push_back(BinomialPP(i, c));
        }
        if (mod > 1) dp.push_back(BinomialPP(mod, 1));
    }

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
        inv = (x + m) % m;
        return inv;
    }

    uint64_t chinese_remainder(const vector<uint32_t>& rem, const vector<uint32_t>& mods){
        uint64_t x, y, res = 0, M = 1;
        for (uint32_t i = 0; i < rem.size(); i++) M *= mods[i];
        for (uint32_t i = 0; i < rem.size(); i++){
            x = M / mods[i];
            y = mod_inverse(x, mods[i]);
            res = (res + ((x * rem[i] % M) * y)) % M;
        }

        return res;
    }

    int64_t binomial(int64_t n, int64_t k){
        if (k > n || n < 0 || k < 0 || !dp.size()) return 0;
        if (n == k || k == 0) return 1;

        vector<uint32_t> rem, mods;
        for (uint32_t i = 0; i < dp.size(); i++){
            mods.push_back(dp[i].mod);
            rem.push_back(dp[i].binomial(n, k));
        }
        return chinese_remainder(rem, mods);
    }
};

int main(){
    auto start = clock();

    int64_t dp[30][30];
    for (int mod = 1; mod < 30; mod++){
        memset(dp, 0, sizeof(dp));
        auto bin = Binomial(mod);

        for (int n = 0; n < 30; n++){
            for (int k = 0; k <= n; k++){
                if (k == 0 || n == k) dp[n][k] = 1 % mod;
                else dp[n][k] = (dp[n - 1][k - 1] + dp[n - 1][k]) % mod;
            }
        }

        for (int n = 0; n < 30; n++){
            for (int k = 0; k < 30; k++){
                assert(dp[n][k] == bin.binomial(n, k));
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

    fprintf(stderr, "\nTime taken = %0.3f\n", (clock()-start) / (double)CLOCKS_PER_SEC); /// 0.022 s
    return 0;
}
