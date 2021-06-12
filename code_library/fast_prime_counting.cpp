#include <bits/stdc++.h>

using namespace std;

/// Magic constants, optimized to answer prime counting queries for n=10^12 but can be tweaked

const int MAXP = 8;
const int MAXN = 52;
const int MAXM = 512512; /// At least the product of the first MAXP - 1 primes

const int MAX = 100000010;

vector<int> primes;
bitset<MAX> is_prime;
unsigned int prod[MAXP], pi[MAX], dp[MAXN][MAXM];

void sieve(){
    is_prime[2] = true;
    for (int i = 3; i < MAX; i += 2) is_prime[i] = true;

    for (int i = 3; i * i < MAX; i += 2){
        for (int j = i * i; is_prime[i] && j < MAX; j += (i << 1)){
            is_prime[j] = false;
        }
    }

    primes.push_back(-1);
    for (int i = 1; i < MAX; i++){
        pi[i] = pi[i - 1] + is_prime[i];
        if (is_prime[i]) primes.push_back(i);
    }
}

void gen(){
    int i, j;
    assert(MAXN >= MAXP);

    sieve();
    for (prod[0] = 1, i = 1; i < MAXP; i++){
        prod[i] = prod[i - 1] * primes[i];
    }

    for (i = 0; i < MAXM; i++) dp[0][i] = i;
    for (i = 1; i < MAXN; i++){
        for (j = 1; j < MAXM; j++){
            dp[i][j] = dp[i - 1][j] - dp[i - 1][j / primes[i]];
        }
    }
}

unsigned long long phi(long long m, int n){
    if (!n) return m;
    if (n < MAXN && m < MAXM) return dp[n][m];
    if (n < MAXP) return dp[n][m % prod[n]] + (m / prod[n]) * dp[n][prod[n]];
    if (primes[n] * primes[n] >= m && m < MAX) return pi[m] - n + 1;

    if ((long long)primes[n] * primes[n] * primes[n] >= m){
        int lim = pi[(int)sqrt(0.5 + m)];
        unsigned long long res = pi[m] - (lim + n - 2) * (lim - n + 1) / 2;
        for (int i = n + 1; i <= lim; i++) res += pi[m / primes[i]];
        return res;
    }

    return phi(m, n - 1) - phi(m / primes[n], n - 1);
}

/***
 *
 * Prime counting function in sublinear time with Meissel-Lehmer algorithm
 * Returns the number of primes not exceeding n
 * Complexity: ~O(n^(2/3))
 *
***/

unsigned long long lehmer(long long n){
    if (n < MAX) return pi[n];

    int s = sqrt(0.5 + n);
    int a = pi[(int)sqrt(0.5 + s)], b = pi[s], c = pi[(int)cbrt(0.5 + n)];
    unsigned long long res = phi(n, a) + (unsigned long long)(b + a - 2) * (b - a + 1) / 2;

    for (int i = a + 1; i <= b; i++){
        unsigned long long w = n / primes[i];
        res -= lehmer(w);
        if (i <= c){
            for (int j = pi[(int)sqrt(0.5 + w)]; j >= i; j--){
                res = res + j - pi[w / primes[j]] - 1;
            }
        }
    }

    return res;
}

int main(){
    auto start = clock();
    gen();
    fprintf(stderr, "Pre-process time = %0.6f\n\n", (clock()-start) / (double)CLOCKS_PER_SEC);  /// 0.295902

    start = clock();

    assert(lehmer(7) == 4);
    assert(lehmer(1000) == 168);
    assert(lehmer(1000000) == 78498);
    assert(lehmer(1000000000) == 50847534);
    assert(lehmer(1e12) == 37607912018LL);
    assert(lehmer(1e13) == 346065536839LL);

    fprintf(stderr, "\nCalculation time = %0.6f\n", (clock()-start) / (double)CLOCKS_PER_SEC);  /// 0.986909
    return 0;
}
