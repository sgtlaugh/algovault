/***
 *
 * Optimized Sieve of Eratosthenes
 *
 * Sieves all numbers from 1 to MAX and stores all the primes in primes[] array
 * The small primes are generated first using a simpler sieve
 * The numbers are chunked into blocks of fixed sizes
 * Each block is processed separately to make it cache-friendly
 * The is_composite[] array is a compressed bit-vector denoting the numbers crossed out for each block
 * The sieve uses a wheel of size 15015 (3*5*7*11*13) to process each block efficiently
 *
 * The algorithm can generate all the prime numbers from 1 to 2^31 in a little under 1 seconds in a 4.00GHz core-i7 PC when compiled with -O2
 * Runtime in CodeForces - 1500 ms with GNU G++ 17
 *
***/

#include <bits/stdc++.h>

#define MAX 2147483647

using namespace std;

const uint32_t block_size = 1048576;

uint32_t s, prime_cnt, sq[65536], sp[65536], primes[106000000];
uint64_t wheel[15015], is_composite[8192], mask[12][62][8192];

inline void setbit(uint64_t* ar, uint32_t bit){
    ar[bit >> 6] |= (1ULL << (bit & 63));
}

inline uint32_t get_idx(uint32_t i, uint32_t j){
    if (sq[j] > i) return (sq[j] - i) >> 1;
    uint32_t x = sp[j] - i % sp[j];
    if ((x & 1) ^ 1) x += sp[j];
    return x >> 1;
}

void small_sieve(){
    for (uint32_t i = 2; i * i < 65536; i++){
        for (uint32_t j = i * i; j < 65536 && !sp[i]; j += i){
            sp[j] = 1;
        }
    }
    for (uint32_t i = 2; i < 65536; i++){
        if (!sp[i]) sp[s] = i, sq[s++] = i * i;
    }
}

void process_block(uint32_t i){
    uint32_t j, k, l, d, m, x, lim = i + block_size, idx = i % 15015, chunk = 0;

    idx = (idx + ((idx * 105) & 127) * 15015) >> 7;
    for (j = 0; (j << 7) < block_size; j += chunk, idx = 0){
        chunk = min(15015 - idx, (block_size >> 7) - j);
        memcpy(is_composite + j, wheel + idx, sizeof(uint64_t) * chunk);
    }
    if (!i) is_composite[0] = (is_composite[0] | 1) & ~110;

    l = block_size >> 1, m = block_size >> 7;
    for (j = 6; j < 18 && i; j++){
        for (x = get_idx(i, j), k = 0, d = j - 6; k < m; k++){
            is_composite[k] |= mask[d][x][k];
        }
    }

    for (j = (i == 0) ? 6 : 18; j < s && sq[j] < lim; j++){
        for (x = get_idx(i, j); x < l; x += sp[j]){
            setbit(is_composite, x);
        }
    }
}

void populate_primes(uint32_t i){
    for (uint32_t j = 0; (j << 7) < block_size; j++){
        uint64_t x = ~is_composite[j];
        while (x){
            uint32_t p = i + (j << 7) + (__builtin_ctzll(x) << 1) + 1;
            if (p <= MAX) primes[prime_cnt++] = p;
            x ^= (-x & x);
        }
    }
}

void fast_sieve(){
    small_sieve();

    for (uint32_t i = 1; i <= 5; i++){
        for (uint32_t j = i + (i > 3); j < 960960; j += sp[i]){
            setbit(wheel, j);
        }
    }

    for (uint32_t i = 6; i <= 17; i++){
        for (uint32_t j = 0; j < sp[i]; j++){
            for (uint32_t k = j; k < (block_size >> 1); k += sp[i]){
                setbit(mask[i - 6][j], k);
            }
        }
    }

    primes[prime_cnt++] = 2;
    for (uint32_t i = 0; i <= MAX; i += block_size){
        process_block(i);
        populate_primes(i);
    }
}

int main(){
    auto start = clock();

    fast_sieve();
    assert(prime_cnt == 105097565);

    vector<int> first_5_primes, last_5_primes;
    for (int i = 0; i < 5; i++){
        first_5_primes.push_back(primes[i]);
        last_5_primes.push_back(primes[prime_cnt - i - 1]);
    }

    assert(first_5_primes == vector<int>({2, 3, 5, 7, 11}));
    assert(last_5_primes == vector<int>({2147483647, 2147483629, 2147483587, 2147483579, 2147483563}));

    fprintf(stderr, "\nTime taken = %0.3f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));   /// Time taken = 0.952
    return 0;
}
