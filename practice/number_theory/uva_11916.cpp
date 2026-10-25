/***
 * Generalized discrete logarithm with Shank's Baby step giant step algorithm
 *
 * discrete_log(a, b, mod): any int a and b (negative ones too), 1 <= mod <= INT_MAX
 * Returns the smallest non-negative x such that a^x = b (modulo mod)
 * If no solution exists, returns -1
 * Works even if a, b and mod are not pairwise co-primes
 * a^0 = 1 (including 0^0), so x = 0 whenever b = 1 (modulo mod), which is every b when mod = 1
 *
 * Complexity: O(sqrt(mod)) expected
 *
***/

#include <bits/stdtr1c++.h>
#include <ext/pb_ds/assoc_container.hpp>

using namespace std;

/// gp_hash_table masks the low bits, and residues like b * 3^i mod 2^30 all share them
struct mix_hash{
    size_t operator()(int x) const{
        return (unsigned long long)x * 0x9E3779B97F4A7C15ULL >> 32;
    }
};

int discrete_log(int a, int b, int mod){
    __gnu_pbds::gp_hash_table<int, int, mix_hash> mp;
    int i, v, x, e = 1, n = sqrt(mod + 0.5) + 1;
    auto expo = [mod](long long base, int k){
        long long res = 1 % mod;
        for (; k; k >>= 1, base = base * base % mod){
            if (k & 1) res = res * base % mod;
        }
        return (int)res;
    };

    a %= mod, b %= mod;
    if (a < 0) a += mod;
    if (b < 0) b += mod;
    if (b == 1 % mod) return 0;

    for (i = 0; i < (n + 3); i++){
        if (e == b) return i;
        mp[(long long)b * e % mod] = i;
        e = (long long)e * a % mod;
    }

    v = e = expo(a, n);
    for (i = 2; i < (n + 3); i++){
        e = (long long)e * v % mod;
        auto it = mp.find(e);
        if (it != mp.end()){
            x = (((long long)n * i - it->second) % mod + mod) % mod;
            return (expo(a, x) == b) ? x : -1;
        }
    }

    return -1;
}

const int MOD = 100000007;

long long pow_mod(long long x, long long e){
    long long res = 1;
    for (x %= MOD; e; e >>= 1, x = x * x % MOD){
        if (e & 1) res = res * x % MOD;
    }
    return res;
}

int main(){
    int t, n, k, b, r;

    scanf("%d", &t);
    for (int cs = 1; cs <= t; cs++){
        scanf("%d %d %d %d", &n, &k, &b, &r);
        set<pair<int, int>> blocked;
        int max_row = 1;
        for (int i = 0, x, y; i < b; i++){
            scanf("%d %d", &x, &y);
            blocked.insert({x, y});
            max_row = max(max_row, x);
        }

        auto ways = [&](long long m){
            long long below = (m - 1) * n, top = n;
            for (auto e : blocked){
                if (e.first == 1) top--;
                else below--;
                if (e.first + 1 <= m && !blocked.count({e.first + 1, e.second})) below--, top++;
            }
            return pow_mod(k - 1, below) * pow_mod(k, top) % MOD;
        };

        int res = max_row;
        if (ways(max_row) != r){
            long long base = ways(max_row + 1), rest = r * pow_mod(base, MOD - 2) % MOD;
            res = max_row + 1 + discrete_log(pow_mod(k - 1, n), rest, MOD);
        }
        printf("Case %d: %d\n", cs, res);
    }
    return 0;
}
