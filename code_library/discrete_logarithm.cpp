/***
 * Generalized discrete logarithm with Shank's Baby step giant step algorithm
 *
 * Given three integers a, b and mod
 * Returns the smallest non-negative x such that (g^x) % p = h % p
 * If no solution exists, returns -1
 * Works even if a, b and mod are not pairwise co-primes
 *
 * Complexity: O(sqrt(mod))
 * Use a custom hash table for more speed
 *
 * Note that it is assumed 0^0 is undefined although Python treats 0^0 as 1 (https://en.wikipedia.org/wiki/Zero_to_the_power_of_zero)
 * If it is required otherwise, should be handled explicitly (Just adding `if b == 1 or mod == 1: return 0` should suffice)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

int expo(long long x, int n, int mod){
    x %= mod;
    long long res = 1;

    while (n){
        if (n & 1) res = res * x % mod;
        x = x * x % mod;
        n >>= 1;
    }

    return res % mod;
}

int discrete_log(int a, int b, int mod){
    tr1::unordered_map<int, int> mp;
    int i, v, x, e = 1, n = sqrt(mod + 0.5) + 1;

    for (i = 0; i < (n + 3); i++){
        if (e == b) return i;
        mp[(long long)b * e % mod] = i;
        e = (long long)e * a % mod;
    }

    v = e = expo(a, n, mod);
    for (i = 2; i < (n + 3); i++){
        e = (long long)e * v % mod;
        if (mp.count(e)){
            x = (((long long)n * i - mp[e]) % mod + mod) % mod;
            return (expo(a, x, mod) == b) ? x : -1;
        }
    }

    return -1;
}

int main(){
    assert(discrete_log(0, 0, 1) == 1);
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

    return 0;
}
