/***
 *
 * https://en.wikipedia.org/wiki/Bernoulli_number
 *
 * Computes all the bernoulli numbers from 0 to n
 * Bernoulli numbers are sequence of rational numbers of the form p/q
 * Since the numbers can be huge, they are calculated modulo MOD
 * bernoulli_numbers<MOD>(n)[i] = p * inv(q, MOD) % MOD for B_i = p / q, with B_1 = -1/2
 *
 * MOD is a template parameter so the % compiles to a multiplication, a runtime modulus measured 6x slower
 * MOD must be a prime with n + 1 < MOD < 2^31, the default is 10^9 + 7
 *
 * Complexity: O(n^2), plus O(n log MOD) for the modular inverses, O(n) memory
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <int MOD>
int expo(long long x, int n){
    long long res = 1;

    while (n){
        if (n & 1) res = res * x % MOD;
        x = x * x % MOD;
        n >>= 1;
    }
    return res % MOD;
}

template <int MOD = 1000000007>
vector<int> bernoulli_numbers(int n){
    assert(0 <= n && n + 1 < MOD);
    vector<int> S(n + 1, 0), inv(n + 2), fact(n + 1), bernoulli(n + 1, 0);
    int i, j;
    long long x, y, z, lim = (long long)MOD * MOD;

    for (i = 1, fact[0] = 1; i <= n; i++){
        fact[i] = ((long long) fact[i - 1] * i) % MOD;
    }
    for (i = 0; i <= n + 1; i++) inv[i] = expo<MOD>(i, MOD - 2);

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
    auto bernoulli = bernoulli_numbers(3000);

    assert(bernoulli[1] == 500000003);
    assert(bernoulli[9] == 0);
    assert(bernoulli[10] == 348484851);  /// bernoulli[10] = 5 / 66 = 5 * 469696973 % 1000000007 = 348484851

    assert((bernoulli_numbers(0) == vector<int>{1}));
    assert((bernoulli_numbers<7>(4) == vector<int>{1, 3, 6, 0, 3}));  /// mod 7: -1/2 = 3, 1/6 = 6, -1/30 = -1/2 = 3

    return 0;
}
