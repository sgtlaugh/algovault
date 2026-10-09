/***
 *
 * Integer Root
 * Exact floor of the square root, cube root and k-th root of 64-bit integers
 *
 * Complexity: O(1) for isqrt and icbrt, O(log n) for iroot
 *
 * isqrt(n): largest r with r * r <= n, for any unsigned long long n
 * icbrt(n): largest r with r * r * r <= n, for any unsigned long long n
 * iroot(n, k): largest r with r^k <= n, for k >= 1
 *
 * sqrtl / cbrtl alone can be off by one near the top of the range, e.g. (long long)sqrtl(1e18 - 1) may give 1e9
 * The floating estimate here is only a starting point that is then corrected exactly
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// True if r^k <= n, without overflow
bool power_at_most(unsigned long long r, int k, unsigned long long n){
    unsigned long long res = 1;
    for (int i = 0; i < k; i++){
        if (r && res > n / r) return false;
        res *= r;
    }
    return res <= n;
}

unsigned long long iroot(unsigned long long n, int k){
    assert(k >= 1);
    if (k == 1 || n < 2) return n;
    unsigned long long r = powl((long double)n, 1.0L / k);
    while (r && !power_at_most(r, k, n)) r--;
    while (power_at_most(r + 1, k, n)) r++;
    return r;
}

unsigned long long isqrt(unsigned long long n){
    return iroot(n, 2);
}

unsigned long long icbrt(unsigned long long n){
    return iroot(n, 3);
}

int main(){
    assert(isqrt(0) == 0 && isqrt(1) == 1 && isqrt(3) == 1 && isqrt(4) == 2 && isqrt(99) == 9 && isqrt(100) == 10);
    assert(isqrt(999999999999999999ULL) == 999999999ULL);
    assert(isqrt(1000000000000000000ULL) == 1000000000ULL);
    assert(isqrt(ULLONG_MAX) == 4294967295ULL);
    assert(isqrt(4294967295ULL * 4294967295ULL) == 4294967295ULL);
    assert(isqrt(4294967295ULL * 4294967295ULL - 1) == 4294967294ULL);

    assert(icbrt(0) == 0 && icbrt(7) == 1 && icbrt(8) == 2 && icbrt(26) == 2 && icbrt(27) == 3);
    assert(icbrt(999999999999999999ULL) == 999999ULL && icbrt(1000000000000000000ULL) == 1000000ULL);
    assert(icbrt(ULLONG_MAX) == 2642245ULL);

    assert(iroot(1024, 10) == 2 && iroot(1023, 10) == 1 && iroot(ULLONG_MAX, 64) == 1 && iroot(ULLONG_MAX, 63) == 2);
    assert(iroot(12345, 1) == 12345 && iroot(1, 50) == 1 && iroot(0, 7) == 0);
    return 0;
}
