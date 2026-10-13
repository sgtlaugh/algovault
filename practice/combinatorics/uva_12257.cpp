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

const int MOD = 1000000007;

Combinatorics comb(1000, MOD);
vector<int> adj[1001];
int subtree[1001];

/// A supervisor stands first in his subtree's block, the children's blocks interleave freely in the remaining places
long long ways(int u){
    long long res = 1;
    subtree[u] = 1;
    for (int v: adj[u]){
        res = res * ways(v) % MOD;
        subtree[u] += subtree[v];
    }

    int free_places = subtree[u] - 1;
    for (int v: adj[u]){
        res = res * comb.nCr(free_places, subtree[v]) % MOD;
        free_places -= subtree[v];
    }
    return res;
}

int main(){
    int t;
    scanf("%d", &t);
    for (int tc = 1; tc <= t; tc++){
        int n;
        scanf("%d", &n);

        vector<bool> has_boss(n + 1, false);
        for (int i = 1; i <= n; i++) adj[i].clear();
        for (int i = 1; i < n; i++){
            int a, b;
            scanf("%d %d", &a, &b);
            adj[a].push_back(b);
            has_boss[b] = true;
        }

        int root = 1;
        while (has_boss[root]) root++;
        printf("Case %d: %lld\n", tc, ways(root));
    }
    return 0;
}
