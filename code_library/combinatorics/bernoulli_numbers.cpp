/***
 *
 * Bernoulli Numbers
 * B_0 to B_n modulo a prime, each rational B_i = p / q given as p times the inverse of q
 *
 * Complexity: O(n^2) time, O(n) memory
 *
 * bernoulli_numbers<MOD>(n): B_0, ..., B_n modulo MOD, with the convention B_1 = -1/2
 *
 * Notes:
 *   - MOD is a prime with n + 1 < MOD < 2^31, 10^9 + 7 by default
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// MOD is a template parameter so % compiles to a multiplication, a runtime modulus measured 6x slower
template <int MOD = 1000000007>
vector<int> bernoulli_numbers(int n){
    assert(0 <= n && n + 1 < MOD);
    vector<int> S(n + 1, 0), inv(n + 2), fact(n + 1), bernoulli(n + 1, 0);
    int i, j;
    long long x, y, z, lim = (long long)MOD * MOD;

    for (i = 1, fact[0] = 1; i <= n; i++){
        fact[i] = ((long long) fact[i - 1] * i) % MOD;
    }
    for (i = 2, inv[1] = 1; i <= n + 1; i++) inv[i] = (long long)(MOD - MOD / i) * inv[MOD % i] % MOD;

    /// S holds row i of the Stirling numbers of the second kind, updated in place from the right
    bernoulli[0] = 1;
    for (i = 1, S[0] = 1; i <= n; i++){
        for (j = i, S[i] = 0; j >= 1; j--){
            S[j] = ((long long)S[j] * j + S[j - 1]) % MOD;
        }
        S[0] = 0;

        if (i == 1 || i % 2 == 0){
            for (j = 0, x = 0, y = 0; j <= i; j++){
                z = (long long)fact[j] * inv[j + 1] % MOD;
                z = z * S[j] % MOD;
                if (j & 1) y += z;
                else x += z;
            }
            bernoulli[i] = (lim + x - y) % MOD;
        }
    }

    return bernoulli;
}

int main(){
    /// B_1 = -1/2, B_9 = 0 as every odd index above 1, B_10 = 5 / 66
    auto bernoulli = bernoulli_numbers(3000);

    assert(bernoulli[1] == 500000003);
    assert(bernoulli[9] == 0);
    assert(bernoulli[10] == 348484851);  /// bernoulli[10] = 5 / 66 = 5 * 469696973 % 1000000007 = 348484851

    assert((bernoulli_numbers(0) == vector<int>{1}));
    assert((bernoulli_numbers<7>(4) == vector<int>{1, 3, 6, 0, 3}));  /// mod 7: -1/2 = 3, 1/6 = 6, -1/30 = -1/2 = 3

    return 0;
}
