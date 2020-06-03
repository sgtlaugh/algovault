/***
 *
 * Binomial co-efficients (n choose k) modulo m
 * Let m = p1^q1 * p2^q2 * ... * pk^qk, then pi^qi must not exceed MAXP for all i in [1, k]
 * Otherwise the function can sometimes calculate the answer if it's 0 and at other times will raise an assertion error
 * This is a generalization of the Lucas theorem, which calculates binomial mod p^q
 * Can pre-process dp[i] for queries if m is constant
 *
***/

#include <stdio.h>
#include <bits/stdtr1c++.h>

#define MAXP 10000010

using namespace std;

namespace crt{
    long long extended_gcd(long long a, long long b, long long& x, long long& y){
        if (!b){
            y = 0, x = 1;
            return a;
        }

        long long g = extended_gcd(b, a % b, y, x);
        y -= ((a / b) * x);
        return g;
    }

    long long mod_inverse(long long a, long long m){
        long long x, y, inv;
        extended_gcd(a, m, x, y);
        inv = (x + m) % m;
        return inv;
    }

    long long chinese_remainder(const vector <long long>& ar, const vector <long long>& mods){
        long long x, y, res = 0, M = 1;
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
    long long mod = 0;

    long long trailing(long long x, long long p){
        long long res = 0;
        while (x){
            x /= p;
            res += x;
        }
        return res;
    }

    long long expo(long long x, long long n, long long m){
        x %= m;
        long long res = 1;

        while (n){
            if (n & 1) res = res * x % m;
            x = x * x % m;
            n >>= 1;
        }

        return res % m;
    }

    long long factorial(long long x, long long p){
        long long res = expo(dp[mod - 1], x / mod, mod);
        if (x >= p) res = res * factorial(x / p, p) % mod;
        return res * dp[x % mod] % mod;
    }

    /// binomial co-efficients modulo p^q (p must be a prime and p^q less than MAXP)
    long long binomial(long long n, long long k, long long p, long long q){
        if (k > n) return 0;
        if (n == k || k == 0) return 1;

        int i;
        for (i = 0, mod = 1; i < q; i++) mod *= p;
        long long t = trailing(n, p) - trailing(k, p) - trailing(n - k, p);
        if (t >= q) return 0;

        if (mod >= MAXP){
            throw std::invalid_argument("p^q is not less than MAXP");
        }

        for (dp[0] = 1, i = 1; i < mod; i++){
            dp[i] = (long long)dp[i - 1] * (i % p ? i : 1) % mod;
        }

        long long m = factorial(k, p) * factorial(n - k, p) % mod;
        long long v = factorial(n, p) * expo(m, (mod / p) * (p - 1) - 1, mod) % mod;

        return v * expo(p, t, mod) % mod;
    }

    long long binomial(long long n, long long k, long long m){
        if (k > n || m == 1) return 0;
        if (n == k || k == 0) return 1;

        /// can use fast factorization below if required
        vector <pair<int, int>> factors;
        for (long long i = 2; i * i <= m; i = i + 1 + (i & 1)){
            int c = 0;
            while (m % i == 0){
                c++;
                m /= i;
            }
            if (c) factors.push_back(make_pair(i, c));
        }
        if (m > 1) factors.push_back(make_pair(m, 1));

        vector <long long> ar, mods;
        for (int i = 0; i < (int)factors.size(); i++){
            long long x = 1;
            for (int j = 0; j < factors[i].second; j++) x *= factors[i].first;
            mods.push_back(x), ar.push_back(binomial(n, k, factors[i].first, factors[i].second));
        }
        return crt::chinese_remainder(ar, mods);
    }
}

int main(){
    using namespace bin;

    int dp[30][30];
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

    typedef tuple <long long, long long, long long, long long> T;
    const vector <T> test_data = {
        T(10, 5, 253, 252),
        T(123456789012345LL, 1000000009, 997 * 1000003, 140000420),
        T(1000000000000000000LL, 1000000000, 1000003, 0),
        T(1000000000000000000LL, 1000000000, 666666667, 0),
        T(100000000000000000LL, 10000000009LL, 1000000007 * 997LL, 353000002471LL),
        T(1000000000000000000LL, 10000000009LL, 1000000007 * 997LL, 361000002527LL),
        T(1000000000000000000LL, 1000, 1000000007 * 997LL, 625000004375LL),
    };

    bool raised_exception = false;

    try{
        cout << binomial(1000000000000000000LL, 1000, 1000000007 * 666666667LL) << endl;
    }
    catch(std::invalid_argument& e){
        raised_exception = true;
    }

    assert(raised_exception);
    return 0;
}
