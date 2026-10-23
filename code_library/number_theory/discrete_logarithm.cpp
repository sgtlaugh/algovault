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

int main(){
    assert(discrete_log(0, 0, 1) == 0);
    assert(discrete_log(1, 1, 1) == 0);
    assert(discrete_log(2, 1, 3) == 0);
    assert(discrete_log(2, 3, 3) == -1);
    assert(discrete_log(2, 3, 4) == -1);
    assert(discrete_log(2, 3, 6) == -1);
    assert(discrete_log(2, 0, 4) == 2);
    assert(discrete_log(6, 0, 8) == 3);
    assert(discrete_log(6, 8, 16) == 3);
    assert(discrete_log(2, 6, 10) == 4);
    assert(discrete_log(5, 33, 58) == 9);
    assert(discrete_log(3589, 58, 97) == -1);
    assert(discrete_log(3589, 1, 97) == 0);
    assert(discrete_log(1000000000, 666666667, 1000000007) == 942190576);
    assert(discrete_log(-2, 3, 5) == 1);
    assert(discrete_log(2, -1, 5) == 2);
    assert(discrete_log(-1, -1, 7) == 1);
    assert(discrete_log(-3, 1, 1) == 0);
    assert(discrete_log(INT_MIN, INT_MIN, INT_MAX) == 1);

    return 0;
}
