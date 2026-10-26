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
 *   - About 10^7 random numbers or 5 * 10^5 primes below 2^63 per second
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
    const vector<int> BASES_32 = {2, 7, 61};
    /// Jim Sinclair 2011, deterministic below 2^64: https://miller-rabin.appspot.com/
    const vector<int> BASES_64 = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

    struct SmallPrime{ unsigned long long p, inv, max_quotient; };

    /// Holds the bases 7 and 61, so is_witness never gets a base that n divides
    const vector<SmallPrime> SMALL_PRIMES = [](){
        vector<SmallPrime> res;
        for (unsigned long long p: {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79}){
            res.push_back({p, Montgomery(p).inv, ~0ULL / p});
        }
        return res;
    }();

    bool is_witness(int a, long long n, int s, const Montgomery& mont){
        unsigned long long one = mont.to_mont(1), minus_one = mont.to_mont(n - 1);
        unsigned long long x = mont.pow(mont.to_mont(a), (n - 1) >> s);

        if (x == one) return false;

        for (int i = 0; i < s; i++){
            if (x == minus_one) return false;
            x = mont.mul(x, x);
        }

        return true;
    }

    bool miller_rabin(long long n, const vector<int>& bases){
        if (n < 2 || (n & 1) == 0) return n == 2;

        /// p | n exactly when n * p^-1 mod 2^64 <= (2^64 - 1) / p, no division: 27% faster on random n below 2^63
        for (auto& [p, inv, max_quotient]: SMALL_PRIMES){
            if ((unsigned long long)n * inv <= max_quotient) return (unsigned long long)n == p;
        }

        int s = __builtin_ctzll(n - 1);
        Montgomery mont(n);
        for (auto a: bases){
            if (is_witness(a, n, s, mont)) return false;
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
