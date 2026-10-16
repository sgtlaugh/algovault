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

int main(){
    int h, w, n;
    scanf("%d %d %d", &h, &w, &n);
    vector<pair<int, int>> cells(n);
    for (auto& c: cells) scanf("%d %d", &c.first, &c.second);
    sort(cells.begin(), cells.end());
    cells.push_back(make_pair(h, w));

    Combinatorics comb(h + w, MOD);
    auto paths = [&](int r1, int c1, int r2, int c2){
        return comb.nCr(r2 - r1 + c2 - c1, r2 - r1);
    };

    /// dp[i] counts paths to cell i that touch no wall before it; sorted order puts every possible predecessor first
    vector<long long> dp(n + 1);
    for (int i = 0; i <= n; i++){
        dp[i] = paths(1, 1, cells[i].first, cells[i].second);
        for (int j = 0; j < i; j++){
            if (cells[j].second > cells[i].second) continue;
            dp[i] = (dp[i] - dp[j] * paths(cells[j].first, cells[j].second, cells[i].first, cells[i].second)) % MOD;
        }
        if (dp[i] < 0) dp[i] += MOD;
    }
    printf("%lld\n", dp[n]);
    return 0;
}
