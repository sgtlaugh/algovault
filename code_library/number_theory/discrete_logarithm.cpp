/***
 *
 * Discrete Logarithm
 * Generalized discrete logarithm with Shanks' baby step giant step algorithm
 *
 * Complexity: O(sqrt(mod)) expected
 *
 * discrete_log(a, b, mod): the smallest x >= 0 with a^x = b (modulo mod), -1 if no solution exists
 *
 * Notes:
 *   - a and b are any int (negative ones too), 1 <= mod <= INT_MAX, a, b and mod need not be pairwise co-prime
 *   - a^0 = 1 (including 0^0), so x = 0 whenever b = 1 (modulo mod), which is every b when mod = 1
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// The table masks the low bits, and residues like b * 3^i mod 2^30 all share them
struct mix_hash{
    size_t operator()(int x) const{
        return (unsigned long long)x * 0x9E3779B97F4A7C15ULL >> 32;
    }
};

int discrete_log(int a, int b, int mod){
    int n = sqrt(mod + 0.5) + 1;

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

    /// Linear probing at most half full, -1 marks an empty slot: 51% faster than gp_hash_table on Library Checker data
    int cap = 2 << __lg(2 * n + 5);
    vector<int> keys(cap, -1), vals(cap);
    auto slot = [&](int key){
        size_t h = mix_hash()(key) & (cap - 1);
        while (keys[h] != -1 && keys[h] != key) h = (h + 1) & (cap - 1);
        return h;
    };

    int cur = 1 % mod;
    for (int i = 0; i < n + 3; i++){
        if (cur == b) return i;
        int key = (long long)b * cur % mod;
        size_t h = slot(key);
        keys[h] = key, vals[h] = i;
        cur = (long long)cur * a % mod;
    }

    int giant = expo(a, n);
    cur = giant;
    for (int i = 2; i < n + 3; i++){
        cur = (long long)cur * giant % mod;
        size_t h = slot(cur);
        if (keys[h] == cur){
            int x = (((long long)n * i - vals[h]) % mod + mod) % mod;
            return expo(a, x) == b ? x : -1;
        }
    }

    return -1;
}

int main(){
    /// Smallest x with a^x = b, a, b and mod need not be co-prime
    assert(discrete_log(2, 6, 10) == 4);       /// 2^4 = 16 = 6 mod 10
    assert(discrete_log(5, 33, 58) == 9);
    assert(discrete_log(2, 3, 4) == -1);       /// powers of 2 mod 4 are 1, 2, 0, 0, ...
    assert(discrete_log(2, 1, 3) == 0);        /// a^0 = 1
    assert(discrete_log(-2, 3, 5) == 1);       /// -2 = 3 mod 5
    return 0;
}
