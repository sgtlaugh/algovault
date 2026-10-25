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
    return res;
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

const int MOD = 1000000007;

long long bernoulli_denominator(int j){
    if (j == 0) return 1;
    if (j == 1) return 2;

    long long den = 1;
    for (int p = 2; p <= j + 1; p++){
        bool prime = true;
        for (int d = 2; d * d <= p; d++) prime &= p % d != 0;
        if (prime && j % (p - 1) == 0) den *= p;
    }

    return den;
}

int main(){
    const int K = 20;
    vector<int> bernoulli = bernoulli_numbers<MOD>(K);
    vector<long long> num(K + 1), den(K + 1);

    for (int j = 0; j <= K; j++){
        den[j] = bernoulli_denominator(j);
        long long x = (long long)bernoulli[j] * den[j] % MOD;
        num[j] = x > MOD / 2 ? x - MOD : x;
    }
    num[1] = 1;

    long long binom[K + 2][K + 2] = {};
    for (int i = 0; i <= K + 1; i++){
        binom[i][0] = 1;
        for (int j = 1; j <= i; j++) binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j];
    }

    int t;
    scanf("%d", &t);
    for (int tc = 0; tc < t; tc++){
        int k;
        scanf("%d", &k);
        vector<long long> a(k + 2, 0), b(k + 2, 1);
        long long m = 1;
        for (int j = 0; j <= k; j++){
            long long p = binom[k + 1][j] * num[j], q = den[j] * (k + 1), g = __gcd(abs(p), q);
            a[k + 1 - j] = p / g, b[k + 1 - j] = q / g;
            m = m / __gcd(m, b[k + 1 - j]) * b[k + 1 - j];
        }

        if (tc) puts("");
        printf("%lld", m);
        for (int i = k + 1; i >= 0; i--) printf(" %lld", a[i] * (m / b[i]));
        puts("");
    }
    return 0;
}
