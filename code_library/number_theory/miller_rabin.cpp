/***
 *
 * Fast primality check with Miller Rabin
 *
 * Uses the deterministic variant of Miller Rabin
 * For more details, check https://miller-rabin.appspot.com/
 *
 * Complexity: O(24) + O(7 * log n)
 * Or, O(24) + O(4 * log n) when n ≤ INT_MAX
 *
 * For large random numbers not exceeding 2^63, it can process 4*10^6 numbers in one second
 * For large primes not exceeding 2^63, it can process around 5*10^5 numbers in one second
 * Requires __int128 (64-bit GCC or Clang)
 * Embeds a checked copy of pollard_rho.cpp's Montgomery to stay standalone
 *
 * To gain more speed, check the following resources
 *     i) https://people.ksp.sk/~misof/primes/
 *    ii) https://github.com/wizykowski/miller-rabin/blob/master/sprp64.h
 *
 * One hack if we need to gain more speed can be to use only the base 921211727, particularly for numbers greater than INT_MAX
 * It's not guaranteed to provide correct answer in all cases, but its highly likely there won't be cases against it
 * Picking random numbers from 1 to 2^31 10^9 times resulted in 600 mismatches with 921211727
 * Picking random numbers from 1 to 2^63 10^9 times resulted in 0 mismatches
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/***
 *
 * Montgomery multiplication for a fixed odd modulus n, values are stored as x * 2^64 mod n
 * Needs no division, ~6x faster than __int128 % n for primality checks
 * The modulus must be odd, even moduli give wrong results silently
 *
***/

/// BEGIN COPY montgomery from code_library/number_theory/pollard_rho.cpp
struct Montgomery{
    unsigned long long n, inv, r2;

    Montgomery(unsigned long long n) : n(n), inv(1){
        for (int i = 0; i < 6; i++) inv *= 2 - n * inv;  /// Newton iteration, doubles the correct bits of n^-1 mod 2^64 each step
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
/// END COPY montgomery

namespace prm{
    const vector<int> BASES_32 = {2, 3, 5, 7};
    const vector<int> BASES_64 = {2, 450775, 1795265022, 9780504, 28178, 9375, 325};

    const vector<int> SMALL_PRIMES = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 193, 407521, 299210837};

    bool is_probable_composite(int a, long long n, int s, const Montgomery& mont){
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

        for (auto &&p: SMALL_PRIMES){
            if (n == p) return true;
            if (n % p == 0) return false;
        }

        int s = __builtin_ctzll(n - 1);
        Montgomery mont(n);
        for (auto a: bases){
            if (is_probable_composite(a, n, s, mont)) return false;
        }

        return true;
    }

    bool is_prime(long long n){
        if (n < INT_MAX) return miller_rabin(n, BASES_32);
        return miller_rabin(n, BASES_64);
    }
}

int main(){
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
