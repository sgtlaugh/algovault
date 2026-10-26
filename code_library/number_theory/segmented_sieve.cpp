/***
 *
 * Segmented Sieve
 * Primality of every number in a window [L, R] far from zero
 *
 * Complexity:
 *   - O(sqrt(max_r) log log max_r) to build
 *   - O((R - L + 1) log log R + pi(sqrt(R))) time and O(R - L) memory per window
 *
 * SegmentedSieve sieve(max_r): base primes up to sqrt(max_r), for max_r <= 1e14
 * sieve.window(L, R): res[i] = 1 if L + i is prime, for 0 <= L <= R <= max_r and R - L <= 1e7
 *
 * Notes:
 *   - Build once for many windows: at max_r = 1e14 the build takes about 70 ms, as long as one 1e7 wide window
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct SegmentedSieve{
    static constexpr int block_size = 32768;

    long long max_r;
    vector<int> primes;

    /// sqrtl is exact on integers up to 1e14, so the floor needs no correction
    SegmentedSieve(long long max_r) : max_r(max_r){
        assert(0 <= max_r && max_r <= 100000000000000LL);
        int limit = sqrtl((long double)max_r);
        vector<char> composite(limit + 1, 0);

        for (int p = 2; p <= limit; p++){
            if (composite[p]) continue;
            primes.push_back(p);
            for (long long j = (long long)p * p; j <= limit; j += p) composite[j] = 1;
        }
    }

    /// Primes up to block_size sieve one L1-sized block at a time, 1.6x faster on 1e7 wide windows to 5e8, 1.5x at 1e14
    vector<char> window(long long L, long long R) const{
        assert(0 <= L && L <= R && R <= max_r && R - L <= 10000000);
        vector<char> res(R - L + 1, 1);
        vector<long long> next_multiple;
        long long first_hi = min(R, L + block_size - 1);

        for (long long p : primes){
            if (p * p > R) break;
            long long j = max(p * p, (L + p - 1) / p * p);
            long long stop = p <= block_size ? first_hi : R;
            for (; j <= stop; j += p) res[j - L] = 0;
            /// Skipped for a one-block window, whose pushes measured 12% slower on 3000 wide windows below 1e9
            if (p <= block_size && first_hi < R) next_multiple.push_back(j);
        }

        for (long long lo = L + block_size; lo <= R; lo += block_size){
            long long hi = min(R, lo + block_size - 1);
            for (size_t i = 0; i < next_multiple.size(); i++){
                long long j = next_multiple[i];
                for (; j <= hi; j += primes[i]) res[j - L] = 0;
                next_multiple[i] = j;
            }
        }

        for (long long x = L; x <= min(R, 1LL); x++) res[x - L] = 0;
        return res;
    }
};

int main(){
    /// One sieve answers many windows: the primes up to 30, then those in [10^9, 10^9 + 21]
    SegmentedSieve sieve(1000000021);
    vector<char> small = sieve.window(0, 30);
    assert(count(small.begin(), small.end(), 1) == 10);
    assert(small[2]);
    assert(small[29]);
    assert(!small[1]);

    vector<char> window = sieve.window(1000000000, 1000000021);
    assert(count(window.begin(), window.end(), 1) == 3);
    assert(window[7]);
    assert(window[9]);
    assert(window[21]);

    return 0;
}
