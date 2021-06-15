/***
 *
 * Binomial co-efficients (n choose k) modulo mod
 *
 * Let mod = p1^q1 * p2^q2 * ... * pk^qk
 * Then pi^qi must not exceed MAXP for all i in [1, k]
 *
 * Otherwise the function can sometimes calculate the answer if it's 0
 * At other times, it will raise an assertion error
 * This is a generalization of the Lucas theorem, which calculates binomial mod p^q
 *
 * Optimization notes:
   * Can pre-process dp[i] for queries if mod is constant
   * Can use fast factorization to factorize mod if required
 *
***/

#include <bits/stdc++.h>

#define MAXP 10000010

using namespace std;

namespace crt{
    uint64_t extended_gcd(uint64_t a, uint64_t b, uint64_t& x, uint64_t& y){
        if (!b){
            y = 0, x = 1;
            return a;
        }

        uint64_t g = extended_gcd(b, a % b, y, x);
        y -= ((a / b) * x);
        return g;
    }

    uint64_t mod_inverse(uint64_t a, uint64_t m){
        uint64_t x, y, inv;
        extended_gcd(a, m, x, y);
        inv = (x + m) % m;
        return inv;
    }

    uint64_t chinese_remainder(const vector <uint64_t>& ar, const vector <uint64_t>& mods){
        uint64_t x, y, res = 0, M = 1;
        for (int i = 0; i < (int)ar.size(); i++) M *= mods[i];
        for (int i = 0; i < (int)ar.size(); i++){
            x = M / mods[i];
            y = mod_inverse(x, mods[i]);
            res = (res + ((x * ar[i] % M) * y)) % M;
        }

        return res;
    }
}

namespace bin{
    int dp[MAXP];

    uint64_t trailing(uint64_t x, uint64_t p){
        uint64_t res = 0;
        while (x){
            x /= p;
            res += x;
        }
        return res;
    }

    uint64_t expo(uint64_t x, uint64_t n, uint64_t m){
        x %= m;
        uint64_t res = 1;

        while (n){
            if (n & 1) res = res * x % m;
            x = x * x % m;
            n >>= 1;
        }

        return res % m;
    }

    uint64_t fact(uint64_t x, uint64_t p, const uint64_t& mod){
        uint64_t res = expo(dp[mod - 1], x / mod, mod);
        if (x >= p) res = res * fact(x / p, p, mod) % mod;
        return res * dp[x % mod] % mod;
    }

    uint64_t binomial(uint64_t n, uint64_t k, uint64_t p, uint64_t q){
        if (k > n) return 0;
        if (n == k || k == 0) return 1;

        uint32_t i;
        uint64_t mod = 1;
        for (i = 0; i < q; i++) mod *= p;
        uint64_t t = trailing(n, p) - trailing(k, p) - trailing(n - k, p);
        if (t >= q) return 0;

        if (mod >= MAXP){
            throw std::invalid_argument("p^q is not less than MAXP");
        }

        for (dp[0] = 1, i = 1; i < mod; i++){
            dp[i] = (uint64_t)dp[i - 1] * (i % p ? i : 1) % mod;
        }

        uint64_t m = fact(k, p, mod) * fact(n - k, p, mod) % mod;
        uint64_t v = fact(n, p, mod) * expo(m, (mod / p) * (p - 1) - 1, mod) % mod;

        return v * expo(p, t, mod) % mod;
    }

    uint64_t binomial(uint64_t n, uint64_t k, uint64_t mod){
        if (k > n || mod == 1) return 0;
        if (n == k || k == 0) return 1;

        /// can use fast factorization below if required
        vector <pair<int, int>> factors;
        for (uint64_t i = 2; i * i <= mod; i = i + 1 + (i & 1)){
            int c = 0;
            while (mod % i == 0){
                c++;
                mod /= i;
            }
            if (c) factors.push_back(make_pair(i, c));
        }
        if (mod > 1) factors.push_back(make_pair(mod, 1));

        vector <uint64_t> ar, mods;
        for (int i = 0; i < (int)factors.size(); i++){
            uint64_t x = 1;
            for (int j = 0; j < factors[i].second; j++) x *= factors[i].first;
            mods.push_back(x), ar.push_back(binomial(n, k, factors[i].first, factors[i].second));
        }
        return crt::chinese_remainder(ar, mods);
    }
}

int main(){
    using namespace bin;
    auto start = clock();

    uint64_t dp[30][30];
    for (int m = 1; m < 30; m++){
        memset(dp, 0, sizeof(dp));

        for (int n = 0; n < 30; n++){
            for (int k = 0; k <= n; k++){
                if (k == 0 || n == k) dp[n][k] = 1 % m;
                else dp[n][k] = (dp[n - 1][k - 1] + dp[n - 1][k]) % m;
            }
        }

        for (int n = 0; n < 30; n++){
            for (int k = 0; k < 30; k++){
                assert(dp[n][k] == binomial(n, k, m));
            }
        }
    }

    typedef tuple<long long, long long, long long, long long> T;
    const vector<T> test_cases = {
        T(10, 5, 253, 252),
        T(123456789012345LL, 1000000009, 997 * 1000003, 140000420),
        T(1000000000000000000LL, 1000000000, 1000003, 0),
        T(1000000000000000000LL, 1000000000, 666666667, 0),
        T(100000000000000000LL, 10000000009LL, 1000000007 * 997LL, 353000002471LL),
        T(1000000000000000000LL, 10000000009LL, 1000000007 * 997LL, 361000002527LL),
        T(1000000000000000000LL, 1000, 1000000007 * 997LL, 625000004375LL),
    };

    for (auto data: test_cases){
        uint64_t n = get<0>(data), k = get<1>(data), mod = get<2>(data), res = get<3>(data);
        assert(binomial(n, k, mod) == res);
    }

    bool raised_exception = false;
    try{
        cout << binomial(1000000000000000000LL, 1000, 1000000007 * 666666667LL) << endl;
    }
    catch(std::invalid_argument& e){
        raised_exception = true;
    }

    assert(raised_exception);
    fprintf(stderr, "\nTime taken = %0.3f\n", (clock()-start) / (double)CLOCKS_PER_SEC);
    return 0;
}
