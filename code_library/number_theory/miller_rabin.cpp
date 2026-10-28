/***
 *
 * Miller-Rabin
 * Deterministic primality test for 64-bit integers
 *
 * Complexity: O(log n), 3 Miller-Rabin rounds below 4759123141 and 7 above
 *
 * prm::is_prime(n): true if n is prime, for any long long n
 *
 * Notes:
 *   - About 1.3 * 10^7 random numbers or 7.5 * 10^5 primes below 2^63 per second
 *   - Requires __int128 (64-bit GCC or Clang)
 *   - Faster tests:
 *       - https://people.ksp.sk/~misof/primes/
 *       - https://github.com/wizykowski/miller-rabin/blob/master/sprp64.h
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Montgomery form for a fixed odd modulus n, values are stored as x * 2^64 mod n, 8x faster than __int128 % n here
/// BEGIN SHARED montgomery
struct Montgomery{
    unsigned long long n, inv, r2;

    Montgomery(unsigned long long n) : n(n), inv(1){
        for (int i = 0; i < 6; i++) inv *= 2 - n * inv;  /// each Newton step doubles the correct bits of n^-1 mod 2^64
        r2 = -n % n;
        r2 = (unsigned __int128)r2 * r2 % n;
    }

    unsigned long long reduce(unsigned __int128 x) const{
        unsigned long long q = (unsigned long long)x * inv, m = ((unsigned __int128)q * n) >> 64, h = x >> 64;
        return h >= m ? h - m : h + n - m;
    }

    unsigned long long mul(unsigned long long x, unsigned long long y) const{
        return reduce((unsigned __int128)x * y);
    }

    unsigned long long add(unsigned long long x, unsigned long long y) const{
        return (x += y) >= n ? x - n : x;
    }

    unsigned long long to_mont(unsigned long long x) const{
        return mul(x, r2);
    }

    unsigned long long pow(unsigned long long x, unsigned long long e) const{
        unsigned long long res = to_mont(1);
        for (; e; e >>= 1, x = mul(x, x)){
            if (e & 1) res = mul(res, x);
        }

        return res;
    }
};
/// END SHARED montgomery

/// BEGIN SHARED miller_rabin
namespace prm{
    /// Jaeschke 1993, deterministic below 4759123141
    const array<int, 3> BASES_32 = {2, 7, 61};
    /// Jim Sinclair 2011, deterministic below 2^64: https://miller-rabin.appspot.com/
    const array<int, 7> BASES_64 = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

    struct SmallPrime{ unsigned long long p, inv, max_quotient; };

    /// Holds the bases 7 and 61, so n never divides a base and no base becomes 0 in Montgomery form
    const vector<SmallPrime> SMALL_PRIMES = [](){
        vector<SmallPrime> res;
        for (unsigned long long p: {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79}){
            res.push_back({p, Montgomery(p).inv, ~0ULL / p});
        }

        return res;
    }();

    /// x = a^d for odd d with n - 1 = d 2^s: a proves n composite unless x is 1 or reaches -1 within s - 1 squarings
    bool is_witness(unsigned long long x, int s, const Montgomery& mont){
        unsigned long long one = mont.to_mont(1), minus_one = mont.n - one;

        if (x == one) return false;

        for (int i = 0; i < s; i++){
            if (x == minus_one) return false;
            x = mont.mul(x, x);
        }

        return true;
    }

    template<size_t K>
    bool miller_rabin(long long n, const array<int, K>& bases){
        if (n < 2 || (n & 1) == 0) return n == 2;

        /// p | n exactly when n * p^-1 mod 2^64 <= (2^64 - 1) / p, no division: 27% faster on random n below 2^63
        for (auto& [p, inv, max_quotient]: SMALL_PRIMES){
            if ((unsigned long long)n * inv <= max_quotient) return (unsigned long long)n == p;
        }

        int s = __builtin_ctzll(n - 1);
        unsigned long long d = (n - 1) >> s;
        Montgomery mont(n);

        /// bases[0] is 2 in both lists, which alone rejects almost every composite, so the rest mostly run on primes
        if (is_witness(mont.pow(mont.to_mont(bases[0]), d), s, mont)) return false;

        /// The other bases in lockstep so their independent multiplications overlap: 63-bit primes 29% faster
        const int k = K - 1;  /// constant, so the lanes unroll: 31-bit primes 12% faster than a runtime count
        unsigned long long x[k], res[k];
        for (int i = 0; i < k; i++) x[i] = mont.to_mont(bases[i + 1]), res[i] = mont.to_mont(1);

        for (; d; d >>= 1){
            if (d & 1){
                for (int i = 0; i < k; i++) res[i] = mont.mul(res[i], x[i]);
            }
            for (int i = 0; i < k; i++) x[i] = mont.mul(x[i], x[i]);
        }

        for (int i = 0; i < k; i++){
            if (is_witness(res[i], s, mont)) return false;
        }

        return true;
    }

    bool is_prime(long long n){
        if (n < 4759123141LL) return miller_rabin(n, BASES_32);
        return miller_rabin(n, BASES_64);
    }
}
/// END SHARED miller_rabin

int main(){
    /// Primes, then composites, among them squares and products of two large primes
    using namespace prm;

    assert(is_prime(2));
    assert(is_prime(3));
    assert(is_prime(97));
    assert(is_prime(1000003));
    assert(is_prime(2783117));
    assert(is_prime(1000000007));
    assert(is_prime(2147483647));
    assert(is_prime(143305320273842137LL));
    assert(is_prime(701874195430938151LL));

    assert(!is_prime(1));
    assert(!is_prime(4));
    assert(!is_prime(245)); // 5 * 7 * 7
    assert(!is_prime(7745740235689LL)); // 2783117 * 2783117
    assert(!is_prime(2783117019481819LL)); // 2783117 * 1000000007
    assert(!is_prime(666666673000000003L)); // 1000000009 * 666666667
    assert(!is_prime(33144425233627921LL)); // 100003 * 474119 * 699053
    assert(!is_prime(3533656326712328192LL)); // 822743477 * 4294967296

    return 0;
}
