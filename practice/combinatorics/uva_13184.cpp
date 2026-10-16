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

const int MAXN = 5000;
const long long MOD = 1000000007;

int main(){
    vector<int> divs(MAXN + 1, 0);
    for (int i = 1; i <= MAXN; i++){
        for (int j = 2 * i; j <= MAXN; j += i) divs[j]++;
    }
    int max_divs = *max_element(divs.begin(), divs.end());

    /// K beyond the largest proper divisor count changes nothing, so only K <= max_divs gets a table
    Combinatorics comb(max_divs, MOD);
    vector<vector<long long>> edges(max_divs + 1, vector<long long>(MAXN + 1, 0));
    vector<vector<long long>> ways(max_divs + 1, vector<long long>(MAXN + 1, 1));
    for (int k = 1; k <= max_divs; k++){
        for (int i = 1; i <= MAXN; i++){
            edges[k][i] = (edges[k][i - 1] + min(divs[i], k)) % MOD;
            ways[k][i] = ways[k][i - 1] * comb.nCr(divs[i], min(divs[i], k)) % MOD;
        }
    }

    int t, n, k;
    scanf("%d", &t);
    while (t--){
        scanf("%d %d", &n, &k);
        k = min(k, max_divs);
        printf("%lld %lld\n", edges[k][n], ways[k][n]);
    }
    return 0;
}
