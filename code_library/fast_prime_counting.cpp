#include <stdio.h>
#include <bits/stdtr1c++.h>

using namespace std;

/// Magic constants, optimized to answer prime counting queries for n=10^12 but can be tweaked

const int MAXP = 8;
const int MAXN = 52;
const int MAXM = 512512; /// At least the product of the first MAXP - 1 primes

const int MAX = 100000010;
const int BLOCK_SIZE = 1048576;

unsigned long long wheel[15015], is_composite[8192];
unsigned int s, p, prod[10], sq[32768], sp[32768], primes[6000010], pi[MAX], dp[MAXN][MAXM];

void setbit(unsigned long long* ar, int bit){
    ar[bit >> 6] |= (1LL << bit);
}

int isqrt(long long x){
    return sqrt(x + 0.666);
}

int icbrt(long long x){
    return cbrt(x + 0.666);
}

int get_idx(unsigned int i, unsigned int j){
    if (sq[j] > i) return (sq[j] - i) >> 1;
    int x = sp[j] - i % sp[j];
    if ((x & 1) ^ 1) x += sp[j];
    return x >> 1;
}

void small_sieve(){
    for (int i = 2; i * i < 32768; i++){
        for (int j = i * i; j < 32768 && !sp[i]; j += i){
            sp[j] = 1;
        }
    }
    for (int i = 2; i < 32768; i++){
        if (!sp[i]) sp[s] = i, sq[s++] = i * i;
    }
}

void process_block(int i){
    unsigned int j, l, x, lim = i + BLOCK_SIZE, idx = i % 15015, chunk = 0;

    idx = (idx + ((idx * 105) & 127) * 15015) >> 7;
    for (j = 0; (j << 7) < BLOCK_SIZE; j += chunk, idx = 0){
        chunk = min(15015 - idx, (BLOCK_SIZE >> 7) - j);
        memcpy(is_composite + j, wheel + idx, sizeof(unsigned long long) * chunk);
    }
    if (!i) is_composite[0] = (is_composite[0] | 1) & ~110ULL;

    l = BLOCK_SIZE >> 1;
    for (j = 6; j < s && sq[j] < lim; j++){
        for (x = get_idx(i, j); x < l; x += sp[j]) setbit(is_composite, x);
    }
}

void populate_primes(int i){
    for (int j = 0; (j << 7) < BLOCK_SIZE; j++){
        unsigned long long mask = ~is_composite[j];
        while (mask){
            int num = i + (j << 7) + (__builtin_ctzll(mask) << 1) + 1;
            if (num < MAX) pi[num] = 1;
            mask ^= (-mask & mask);
        }
    }
}

void fast_sieve(){
    small_sieve();
    for (int i = 1; i <= 5; i++){
        for (int j = i + (i > 3); j < 960960; j += sp[i]){
            setbit(wheel, j);
        }
    }

    pi[2] = 1;
    for (int i = 0; i <= MAX; i += BLOCK_SIZE){
        process_block(i);
        populate_primes(i);
    }

    for (int i = 2; i < MAX; i++){
        if (pi[i]) primes[++p] = i;
        pi[i] += pi[i - 1];
    }
}

void generate(){
    int i, j;
    fast_sieve();
    for (i = 0; i < MAXM; i++) dp[0][i] = i;
    for (prod[0] = 1, i = 1; i < MAXP; i++) prod[i] = prod[i - 1] * primes[i];

    assert(MAXN >= MAXP);
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
    if (primes[n] * primes[n] >= m) return pi[m] - n + 1;

    if (primes[n] * primes[n] * primes[n] >= m){
        int lim = pi[isqrt(m)];
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

    int s = isqrt(n);
    int a = pi[isqrt(s)], b = pi[s], c = pi[icbrt(n)];
    unsigned long long res = phi(n, a) + (unsigned long long)(b + a - 2) * (b - a + 1) / 2;

    for (int i = a + 1; i <= b; i++){
        unsigned long long w = n / primes[i];
        res -= lehmer(w);
        if (i <= c){
            for (int j = pi[isqrt(w)]; j >= i; j--){
                res = res + j - pi[w / primes[j]] - 1;
            }
        }
    }

    return res;
}

int main(){
    clock_t start = clock();
    generate();
    fprintf(stderr, "Pre-process time = %0.6f\n\n", (clock()-start) / (1.0*CLOCKS_PER_SEC));  /// 0.295902

    cout << lehmer(7) << endl;     /// 4
    cout << lehmer(1e3) << endl;   /// 168
    cout << lehmer(1e6) << endl;   /// 78498
    cout << lehmer(1e9) << endl;   /// 50847534
    cout << lehmer(1e12) << endl;  /// 37607912018

    start = clock();
    cout << lehmer(1e13) << endl;  /// 346065536839
    fprintf(stderr, "\nCalculation time = %0.6f\n", (clock()-start) / (1.0*CLOCKS_PER_SEC));  /// 0.986909

    return 0;
}
