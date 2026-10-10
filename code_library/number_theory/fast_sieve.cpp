/***
 *
 * Optimized Sieve of Eratosthenes
 *
 * FastSieve sieve(n) sieves all numbers from 1 to n <= 2^31 - 1 and stores the primes in ascending order in sieve.primes
 * The small primes are generated first using a simpler sieve
 * The numbers are chunked into blocks of fixed sizes
 * Each block is processed separately to make it cache-friendly
 * The is_composite[] array is a compressed bit-vector denoting the numbers crossed out for each block
 * The sieve uses a wheel of size 15015 (3*5*7*11*13) to process each block efficiently
 *
 * Complexity: O(n log log n) time, O(pi(n)) memory
 *   - primes takes 4 bytes per prime, 420 MB at n = 2^31 - 1, and the masks for the primes 17 to 61 take 49 MB for any n
 *   - The per block bit vector over odd numbers is 64 KB
 *
 * The algorithm can generate all the prime numbers from 1 to 2^31 in a little under 1 seconds in a 4.00GHz core-i7 PC when compiled with -O2
 * Runtime in CodeForces - 1500 ms with GNU G++ 17
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FastSieve{
    static constexpr uint32_t block_size = 1048576;

    uint32_t n, s = 0;
    vector<uint32_t> primes, sq, sp;
    vector<uint64_t> wheel, is_composite, mask;  /// mask[((t - 6) * 62 + x) * 8192 + k] = word k of the multiples of sp[t] in a block starting at offset x, 6 <= t <= 17

    FastSieve(uint32_t n) : n(n), sq(65536), sp(65536), wheel(15015), is_composite(block_size >> 7), mask(12 * 62 * (block_size >> 7)){
        assert(n <= 2147483647u);
        small_sieve();

        for (uint32_t i = 1; i <= 5; i++){
            for (uint32_t j = i + (i > 3); j < 960960; j += sp[i]){
                setbit(wheel.data(), j);
            }
        }

        for (uint32_t i = 6; i <= 17; i++){
            for (uint32_t j = 0; j < sp[i]; j++){
                for (uint32_t k = j; k < (block_size >> 1); k += sp[i]){
                    setbit(&mask[((i - 6) * 62 + j) * (block_size >> 7)], k);
                }
            }
        }

        if (n < 2) return;
        primes.reserve(1.25506 * n / log(n) + 1);  /// pi(x) < 1.25506 x / ln x for x > 1, Rosser and Schoenfeld
        primes.push_back(2);
        for (uint32_t i = 0; i <= n; i += block_size){
            process_block(i);
            populate_primes(i);
        }
    }

    static void setbit(uint64_t* ar, uint32_t bit){
        ar[bit >> 6] |= (1ULL << (bit & 63));
    }

    uint32_t get_idx(uint32_t i, uint32_t j) const{
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
        uint64_t* bits = is_composite.data();

        idx = (idx + ((idx * 105) & 127) * 15015) >> 7;
        for (j = 0; (j << 7) < block_size; j += chunk, idx = 0){
            chunk = min(15015 - idx, (block_size >> 7) - j);
            memcpy(bits + j, wheel.data() + idx, sizeof(uint64_t) * chunk);
        }
        if (!i) bits[0] = (bits[0] | 1) & ~110;

        l = block_size >> 1, m = block_size >> 7;
        for (j = 6; j < 18 && i; j++){
            const uint64_t* pre = &mask[((j - 6) * 62 + get_idx(i, j)) * m];
            for (k = 0; k < m; k++){
                bits[k] |= pre[k];
            }
        }

        for (j = (i == 0) ? 6 : 18; j < s && sq[j] < lim; j++){
            for (x = get_idx(i, j), d = sp[j]; x < l; x += d){
                setbit(bits, x);
            }
        }
    }

    /// Stops at the word holding n and trims past it, a per prime p <= n check measured 4% slower overall at n = 2^31 - 1
    void populate_primes(uint32_t i){
        const uint64_t* bits = is_composite.data();
        uint32_t words = min<uint64_t>(block_size, (uint64_t)n - i + 128) >> 7;

        for (uint32_t j = 0; j < words; j++){
            uint64_t x = ~bits[j];
            while (x){
                primes.push_back(i + (j << 7) + (__builtin_ctzll(x) << 1) + 1);
                x ^= (-x & x);
            }
        }
        while (primes.back() > n) primes.pop_back();
    }
};

int main(){
    auto start = clock();

    FastSieve sieve(2147483647);
    const vector<uint32_t>& primes = sieve.primes;
    assert(primes.size() == 105097565);

    vector<int> first_5_primes, last_5_primes;
    for (int i = 0; i < 5; i++){
        first_5_primes.push_back(primes[i]);
        last_5_primes.push_back(primes[primes.size() - i - 1]);
    }

    assert(first_5_primes == vector<int>({2, 3, 5, 7, 11}));
    assert(last_5_primes == vector<int>({2147483647, 2147483629, 2147483587, 2147483579, 2147483563}));

    assert(FastSieve(0).primes.empty() && FastSieve(1).primes.empty());
    assert(FastSieve(2).primes == vector<uint32_t>({2}));
    assert(FastSieve(30).primes == vector<uint32_t>({2, 3, 5, 7, 11, 13, 17, 19, 23, 29}));
    assert(FastSieve(1000000).primes.size() == 78498);

    fprintf(stderr, "\nTime taken = %0.3f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));   /// Time taken = 0.952
    return 0;
}
