#include <stdio.h>
#include <bits/stdtr1c++.h>


using namespace std;

typedef __int128 int128;

/// Magic constants, optimized to answer prime counting queries for n=10^12 but can be tweaked

const int MAXN = 50;
const int MAXM = 1000010;

const int MAX = 100000010;
const int BLOCK_SIZE = 1048576;

unsigned int s, p, sq[32768], sp[32768], primes[6000010], pi[MAX];
unsigned long long wheel[15015], is_composite[8192], pi_sum[MAX], dp[MAXN][MAXM];

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
            if (num < MAX) pi[num] = 1, pi_sum[num] = num;
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

    pi[2] = 1, pi_sum[2] = 2;
    for (int i = 0; i <= MAX; i += BLOCK_SIZE){
        process_block(i);
        populate_primes(i);
    }

    for (int i = 2; i < MAX; i++){
        if (pi[i]) primes[++p] = i;
        pi[i] += pi[i - 1], pi_sum[i] += pi_sum[i - 1];
    }
}

void generate(){
    fast_sieve();

    for (int i = 0; i < MAXM; i++) dp[0][i] = (long long)i * (i + 1) / 2;
    for (int i = 1; i < MAXN; i++){
        for (int j = 1; j < MAXM; j++){
            dp[i][j] = dp[i - 1][j] - dp[i - 1][j / primes[i]] * primes[i];
        }
    }
}

int128 phi(long long m, int n){
    if (!n) return (int128)m * (m + 1) / 2;
    if (n < MAXN && m < MAXM) return dp[n][m];
    if ((long long)primes[n] * primes[n] >= m) return pi_sum[m] - pi_sum[primes[n]] + 1;
    return phi(m, n - 1) - phi(m / primes[n], n - 1) * primes[n];
}

/***
 *
 * Prime counting function in sublinear time with Meissel-Lehmer algorithm
 * Returns the sum of all the primes not exceeding n
 * Complexity: ~O(n^(2/3))
 *
***/

int128 lehmer(long long n){
    if (n < MAX) return pi_sum[n];


    unsigned int s = isqrt(n), c = icbrt(n);
    int128 res = phi(n, pi[c]) + pi_sum[c] - 1;

    for (int i = pi[c] + 1; primes[i] <= (s + 1); i++){
        int128 w = lehmer(n / primes[i]) - pi_sum[primes[i] - 1];
        res -= w * primes[i];
    }

    return res;
}

/// Just for demonstration, not handling the case for v=0 intentionally
void print128(int128 v){
    if (!v){
        puts("");
        return;
    }
    print128(v / 10);
    printf("%d", (int)(v % 10));
}

int main(){
    clock_t start = clock();
    generate();
    fprintf(stderr, "Pre-process time = %0.6f\n\n", (clock()-start) / (1.0*CLOCKS_PER_SEC));  /// 0.433484

    cout << (long long)lehmer(7) << endl;     /// 17
    cout << (long long)lehmer(1e3) << endl;   /// 76127
    cout << (long long)lehmer(1e6) << endl;   /// 37550402023
    cout << (long long)lehmer(1e9) << endl;   /// 24739512092254535
    cout << (long long)lehmer(1e10) << endl;  /// 2220822432581729238

    start = clock();
    int128 res = lehmer(1e13);
    print128(res);                            /// 1699246443377779418889494

    fprintf(stderr, "\nCalculation time = %0.6f\n", (clock()-start) / (1.0*CLOCKS_PER_SEC));  /// 1.486966
    return 0;
}
