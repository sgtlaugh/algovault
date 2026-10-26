/***
 *
 * Fast Sieve
 * Optimized Sieve of Eratosthenes for every prime up to 2^31 - 1
 *
 * Complexity: O(n log log n) time, about 550 MB virtual and 460 MB resident memory at n = 2^31 - 1
 *
 * FastSieve sieve(n): sieves all numbers from 1 to n <= 2^31 - 1
 * sieve.primes: the primes up to n in ascending order
 *
 * Notes:
 *   - Every prime below 2^31 in under 1 s on a 4.00 GHz Core i7 at -O2, 1500 ms on Codeforces with GNU G++ 17
 *   - Memory is 4 bytes per prime plus 49 MB of masks once n >= 2^20, under 1 MB more below that
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct FastSieve{
    static constexpr uint32_t block_size = 1048576;

    uint32_t n, small_count = 0;
    vector<uint32_t> primes, squares, small_primes;
    /// mask[((t - 6) * 62 + x) * 8192 + k]: word k of small_primes[t]'s multiples in a block at offset x, 6 <= t <= 17
    vector<uint64_t> wheel, is_composite, mask;

    /// Only blocks after the first read the masks, so n < block_size skips their 49 MB and 34 ms
    FastSieve(uint32_t n) : n(n), squares(65536), small_primes(65536), wheel(15015), is_composite(block_size >> 7), mask(n >= block_size ? 12 * 62 * (block_size >> 7) : 0){
        assert(n <= 2147483647u);
        small_sieve();

        for (uint32_t i = 1; i <= 5; i++){
            for (uint32_t j = i + (i > 3); j < 960960; j += small_primes[i]){
                setbit(wheel.data(), j);
            }
        }

        for (uint32_t i = 6; i <= 17 && !mask.empty(); i++){
            for (uint32_t j = 0; j < small_primes[i]; j++){
                for (uint32_t k = j; k < (block_size >> 1); k += small_primes[i]){
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

    void small_sieve(){
        for (uint32_t i = 2; i * i < 65536; i++){
            for (uint32_t j = i * i; j < 65536 && !small_primes[i]; j += i){
                small_primes[j] = 1;
            }
        }

        for (uint32_t i = 2; i < 65536; i++){
            if (!small_primes[i]) small_primes[small_count] = i, squares[small_count++] = i * i;
        }
    }

    static void setbit(uint64_t* ar, uint32_t bit){
        ar[bit >> 6] |= (1ULL << (bit & 63));
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

        for (j = (i == 0) ? 6 : 18; j < small_count && squares[j] < lim; j++){
            for (x = get_idx(i, j), d = small_primes[j]; x < l; x += d){
                setbit(bits, x);
            }
        }
    }

    uint32_t get_idx(uint32_t i, uint32_t j) const{
        if (squares[j] > i) return (squares[j] - i) >> 1;
        uint32_t x = small_primes[j] - i % small_primes[j];
        if ((x & 1) ^ 1) x += small_primes[j];
        return x >> 1;
    }

    /// Stops at the word holding n and trims past it, a per prime p <= n check measured 4% slower at n = 2^31 - 1
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
    /// The primes up to 30
    FastSieve sieve(30);
    assert(sieve.primes == vector<uint32_t>({2, 3, 5, 7, 11, 13, 17, 19, 23, 29}));

    /// pi(10^6) = 78498, the largest prime below 10^6 is 999983
    FastSieve million(1000000);
    assert(million.primes.size() == 78498);
    assert(million.primes.back() == 999983);

    return 0;
}
