#include "../common.h"

#define main library_main
#include "../../code_library/combinatorics/bernoulli_numbers.cpp"
#undef main

long long power(long long x, long long n, long long mod){
    long long res = 1;
    for (x %= mod; n; n >>= 1, x = x * x % mod) if (n & 1) res = res * x % mod;
    return res;
}

/// The classic recurrence sum C(m + 1, j) B_j = 0, independent of the library's Stirling number formula
vector<int> reference(int n, long long mod){
    vector<vector<long long>> binom(n + 2, vector<long long>(n + 2, 0));
    for (int i = 0; i <= n + 1; i++){
        binom[i][0] = 1;
        for (int j = 1; j <= i; j++) binom[i][j] = (binom[i - 1][j - 1] + binom[i - 1][j]) % mod;
    }

    vector<int> b(n + 1);
    b[0] = 1;
    for (int m = 1; m <= n; m++){
        long long s = 0;
        for (int j = 0; j < m; j++) s = (s + binom[m + 1][j] * b[j]) % mod;
        b[m] = (mod - s) % mod * power(m + 1, mod - 2, mod) % mod;
    }
    return b;
}

int main(){
    const long long MOD = 1000000007;
    const int n = 3008;
    auto expected = reference(n, MOD);
    auto bernoulli = bernoulli_numbers(n);
    assert(bernoulli == expected);

    /// A second call must not see the previous run's values, every shorter n is a prefix
    assert(bernoulli_numbers(n) == expected);
    for (int m = 0; m <= 200; m++) assert(bernoulli_numbers(m) == vector<int>(expected.begin(), expected.begin() + m + 1));

    /// Other primes, up to the largest n the header allows, n + 1 < MOD
    assert(bernoulli_numbers<998244353>(1000) == reference(1000, 998244353));
    assert(bernoulli_numbers<2147483647>(500) == reference(500, 2147483647));
    assert(bernoulli_numbers<1009>(1007) == reference(1007, 1009));
    assert(bernoulli_numbers<3>(1) == reference(1, 3));

    /// Exact small values: B_1 = -1/2, B_2 = 1/6, B_4 = -1/30, B_12 = -691/2730
    assert(bernoulli[1] == MOD - power(2, MOD - 2, MOD));
    assert(bernoulli[2] == power(6, MOD - 2, MOD));
    assert(bernoulli[4] == MOD - power(30, MOD - 2, MOD));
    assert(bernoulli[12] == (MOD - 691) * power(2730, MOD - 2, MOD) % MOD);

    return 0;
}
