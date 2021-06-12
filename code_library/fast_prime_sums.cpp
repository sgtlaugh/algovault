/***
 *
 * Prime counting function in sublinear time with Meissel-Lehmer algorithm
 * Extension to calculate sum of primes using the same idea instead of the count
 *
 * The function prime_sum(n) returns the number of primes not exceeding n
 * It is just a templatized wrapper of lehmer(n)
 *
 * Complexity: Roughly ~O(n^(2/3))
 *
***/

#include <bits/stdc++.h>

using namespace std;

/// Magic constants, optimized to answer prime counting queries for n=10^12 but can be tweaked

const int MAXN = 50;
const int MAXM = 1000010;
const int MAXV = 10000010;

int pi[MAXV];
long long pi_sum[MAXV], dp[MAXN][MAXM];

vector<int> primes;
bitset<MAXV> is_prime;

void sieve(){
    is_prime[2] = true;
    for (int i = 3; i < MAXV; i += 2) is_prime[i] = true;

    for (int i = 3; i * i < MAXV; i += 2){
        for (int j = i * i; is_prime[i] && j < MAXV; j += (i << 1)){
            is_prime[j] = false;
        }
    }

    primes.push_back(-1);
    for (int i = 1; i < MAXV; i++){
        pi[i] = pi[i - 1], pi_sum[i] = pi_sum[i - 1];
        if (is_prime[i]){
            primes.push_back(i);
            pi[i]++, pi_sum[i] += i;
        }
    }
}

void gen(){
    sieve();
    for (int i = 0; i < MAXM; i++) dp[0][i] = (long long)i * (i + 1) / 2;
    for (int i = 1; i < MAXN; i++){
        for (int j = 1; j < MAXM; j++){
            dp[i][j] = dp[i - 1][j] - dp[i - 1][j / primes[i]] * primes[i];
        }
    }
}

typedef uint64_t int128;

template <typename T>
T phi(T m, int n){
    if (!n) return (T)m * (m + 1) / 2;
    if (n < MAXN && m < MAXM) return dp[n][m];
    if (m < MAXV && (uint64_t)primes[n] * primes[n] >= m) return pi_sum[m] - pi_sum[primes[n]] + 1;
    return phi(m, n - 1) - phi(m / primes[n], n - 1) * primes[n];
}

template <typename T>
T lehmer(T n){
    if (n < MAXV) return pi_sum[n];

    int s = sqrt(0.5 + n), c = cbrt(0.5 + n);
    T res = phi(n, pi[c]) + pi_sum[c] - 1;

    for (int i = pi[c] + 1; primes[i] <= (s + 1); i++){
        T w = lehmer(n / primes[i]) - pi_sum[primes[i] - 1];
        res -= w * primes[i];
    }

    return res;
}

__int128 prime_sum(long long n){
    if (n <= UINT_MAX) return lehmer((uint64_t)n);
    return lehmer((__int128)n);
}

int main(){
    auto start = clock();
    gen();
    fprintf(stderr, "Pre-process time = %0.3f\n\n", (clock()-start) / (double)CLOCKS_PER_SEC);  /// 0.399

    start = clock();

    assert(prime_sum(1000) == 76127);
    assert(prime_sum(1000000) == 37550402023);
    assert(prime_sum(1e9) == 24739512092254535LL);
    assert(prime_sum(1e10) == 2220822432581729238LL);
    assert(prime_sum(UINT_MAX) == 425649736193687430LL);

    assert(prime_sum(1e11) == (__int128)603698 * 333721625289043LL);      /// 201467077743744681014
    assert(prime_sum(1e12) == (__int128)15929208151LL * 1157344946327LL); /// 18435588552550705911377
    assert(prime_sum(1e13) == (__int128)10166702 * 167138413556114797LL); /// 1699246443377779418889494

    fprintf(stderr, "\nCalculation time = %0.3f\n", (clock()-start) / (double)CLOCKS_PER_SEC);  /// 3.520
    return 0;
}
