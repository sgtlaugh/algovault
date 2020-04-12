#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define MAX 2000000010
#define min(a,b) ((a)<(b) ? (a):(b))

/***
 *
 * Optimized Sieve of Eratosthenes
 * Sieves all numbers from 1 to MAX and stores all the primes in primes[] array
 * The small primes are generated first using a simpler sieve which are later used in the actual sieve
 * The numbers are chunked into blocks of fixed sizes and each blocks are processed separately to make it cache-friendly
 * The is_composite[] array is a compressed bit-vector which denotes the numbers crossed out for each block
 * The sieve uses a wheel of size 15015 (3*5*7*11*13) and utilizes this to process each block efficiently (https://en.wikipedia.org/wiki/Wheel_factorization)
 * Note the sieve processes all the numbers from 1 to MAX, but might process upto an extra block_size numbers after MAX
 *
 * The algorithm can generate all the prime numbers from 1 to 2*10^9 in a little under 1 seconds in a 4.00GHz core-i7 PC when compiled with -O2
 * Runtime in CF/SPJ - 1600 ms with GNU GCC C11 5.1.0
 *
***/

const int block_size = 1048576;

int s, prime_cnt, sq[65536], sp[65536], primes[100000010];
unsigned long long pow2[64], wheel[15015], is_composite[8192], mask[12][62][8192];

void setbit(unsigned long long* ar, int bit){
    ar[bit >> 6] |= pow2[bit & 63];
}

int get_idx(int i, int j){
    if (sq[j] > i) return (sq[j] - i) >> 1;
    int x = sp[j] - i % sp[j];
    if ((x & 1) ^ 1) x += sp[j];
    return x >> 1;
}

void small_sieve(){
    for (int i = 2; i * i < 65536; i++){
        for (int j = i * i; j < 65536 && !sp[i]; j += i){
            sp[j] = 1;
        }
    }
    for (int i = 2; i < 65536; i++){
        if (!sp[i]) sp[s] = i, sq[s++] = i * i;
    }
}

void process_block(int i){
    int j, k, l, d, m, x, lim = i + block_size, idx = i % 15015, chunk = 0;

    idx = (idx + ((idx * 105) & 127) * 15015) >> 7;
    for (j = 0; (j << 7) < block_size; j += chunk, idx = 0){
        chunk = min(15015 - idx, (block_size >> 7) - j);
        memcpy(is_composite + j, wheel + idx, sizeof(unsigned long long) * chunk);
    }
    if (!i) is_composite[0] = (is_composite[0] | 1) & ~110ULL;

    l = block_size >> 1, m = block_size >> 7;
    for (j = 6; j < 18 && i; j++){
        for (x = get_idx(i, j), k = 0, d = j - 6; k < m; k++) is_composite[k] |= mask[d][x][k];
    }

    for (j = (i == 0) ? 6 : 18; j < s && sq[j] < lim; j++){
        for (x = get_idx(i, j); x < l; x += sp[j]) setbit(is_composite, x);
    }
}

void populate_primes(int i){
    for (int j = 0; (j << 7) < block_size; j++){
        unsigned long long mask = ~is_composite[j];
        while (mask){
            primes[prime_cnt++] = i + (j << 7) + (__builtin_ctzll(mask) << 1) + 1;
            mask ^= (-mask & mask);
        }
    }
}

void fast_sieve(){
    small_sieve();

    for (int i = 0; i < 64; i++) pow2[i] = 1ULL << i;
    for (int i = 1; i <= 5; i++){
        for (int j = i + (i > 3); j < 960960; j += sp[i]){
            setbit(wheel, j);
        }
    }

    for (int i = 6; i <= 17; i++){
        for (int j = 0; j < sp[i]; j++){
            for (int k = j; k < (block_size >> 1); k += sp[i]) setbit(mask[i - 6][j], k);
        }
    }

    primes[prime_cnt++] = 2;
    for (int i = 0; i <= MAX; i += block_size){
        process_block(i);
        populate_primes(i);
    }
}

int main(){
    clock_t start = clock();
    fast_sieve();

    printf("%d\n", prime_cnt);   /// 98254282 (all primes from 1 to 2000683007)
    for (int i = 0; i < 10; i++) printf("%d ", primes[i]); /// 2 3 5 7 11 13 17 19 23 29
    puts("");

    fprintf(stderr, "\nTime taken = %0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));   /// Time taken = 0.990301
    return 0;
}
