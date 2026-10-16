/***
 *
 * Combinatorics
 * Factorial tables, nCr and nPr modulo a prime
 *
 * Complexity: O(n) to build Combinatorics(n, mod), O(1) per query afterwards
 *
 * Combinatorics comb(n, mod): tables for 0..n, mod must be a prime with n < mod < 2^31
 * comb.nCr(n, k), comb.nPr(n, k): 0 when k < 0 or k > n, n must not exceed the table size
 * comb.inverse(i): inverse of i for 1 <= i <= n, comb.factorial(i), comb.inv_factorial(i)
 *
 * Extended gcd, modular inverse for any modulus and diophantine equations: extended_gcd.cpp
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct Combinatorics{
    long long mod;
    vector<long long> fact, inv_fact, inv;

    Combinatorics(int n, long long mod) : mod(mod), fact(n + 1), inv_fact(n + 1), inv(n + 1){
        assert(n >= 0 && n < mod && mod < (1LL << 31));

        fact[0] = 1;
        for (int i = 1; i <= n; i++) fact[i] = fact[i - 1] * i % mod;

        if (n >= 1) inv[1] = 1;
        for (int i = 2; i <= n; i++) inv[i] = mod - (mod / i) * inv[mod % i] % mod;  /// needs a prime mod

        inv_fact[0] = 1;
        for (int i = 1; i <= n; i++) inv_fact[i] = inv_fact[i - 1] * inv[i] % mod;
    }

    long long factorial(int i) const{
        assert(0 <= i && i < (int)fact.size());
        return fact[i];
    }

    long long inv_factorial(int i) const{
        assert(0 <= i && i < (int)fact.size());
        return inv_fact[i];
    }

    long long inverse(int i) const{
        assert(1 <= i && i < (int)inv.size());
        return inv[i];
    }

    long long nCr(int n, int k) const{
        assert(n < (int)fact.size());
        if (k < 0 || k > n) return 0;
        return fact[n] * inv_fact[k] % mod * inv_fact[n - k] % mod;
    }

    long long nPr(int n, int k) const{
        assert(n < (int)fact.size());
        if (k < 0 || k > n) return 0;
        return fact[n] * inv_fact[n - k] % mod;
    }
};

const long long MOD = 1000000007;
const int MAXM = 500000;

int main(){
    Combinatorics comb(MAXM, MOD);
    int p;
    scanf("%d", &p);
    vector<long long> pa(MAXM + 1), pb(MAXM + 1);
    while (p--){
        int m, n;
        scanf("%d %d", &m, &n);
        pa[0] = pb[0] = 1;
        for (int i = 1; i <= m; i++){
            pa[i] = pa[i - 1] * (n - 1) % MOD;
            pb[i] = pb[i - 1] * n % MOD;
        }

        /// the LIS is n exactly when 1..n appears as a subsequence; fix where the greedy match of 1..n ends at i:
        /// the i - n unmatched slots before it dodge the next needed value (n - 1 choices), the rest are free
        long long res = 0;
        for (int i = n; i <= m; i++){
            res = (res + comb.nCr(i - 1, n - 1) * pa[i - n] % MOD * pb[m - i]) % MOD;
        }
        printf("%lld\n", res);
    }
    return 0;
}
