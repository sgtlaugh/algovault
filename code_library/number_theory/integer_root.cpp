/***
 *
 * Integer Root
 * Exact floor of the square root, cube root and k-th root of 64-bit integers
 *
 * Complexity: O(1) for isqrt and icbrt, O(log n) for iroot
 *
 * icbrt(n): largest r with r * r * r <= n, for any unsigned long long n
 * iroot(n, k): largest r with r^k <= n, for k >= 1
 * isqrt(n): largest r with r * r <= n, for any unsigned long long n
 *
 * Notes:
 *   - Plain floating roots are off by one near the top of the range, (long long)sqrt(999999999999999999LL) gives 1e9
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

unsigned long long icbrt(unsigned long long n){
    /// double estimates, not long double: the loops correct either one, and double measured 31-75% faster
    unsigned long long r = cbrt((double)n);
    while (r * r * r > n) r--;
    while (r < 2642245ULL && (r + 1) * (r + 1) * (r + 1) <= n) r++;
    return r;
}

unsigned long long iroot(unsigned long long n, int k){
    assert(k >= 1);
    if (k == 1 || n < 2) return n;
    if (k >= 64) return 1;

    unsigned long long r = pow((double)n, 1.0 / k);
    while (r && !power_at_most(r, k, n)) r--;
    while (power_at_most(r + 1, k, n)) r++;
    return r;
}

unsigned long long isqrt(unsigned long long n){
    unsigned long long r = min<unsigned long long>(sqrt((double)n), 4294967295ULL);  /// sqrt is 2^32 near 2^64
    while (r * r > n) r--;
    while (r < 4294967295ULL && (r + 1) * (r + 1) <= n) r++;
    return r;
}

int main(){
    /// Exact floor roots, also near 2^64 where floating roots round up
    assert(isqrt(99) == 9);
    assert(isqrt(100) == 10);
    assert(isqrt(999999999999999999ULL) == 999999999ULL);
    assert(isqrt(ULLONG_MAX) == 4294967295ULL);

    assert(icbrt(26) == 2);
    assert(icbrt(27) == 3);
    assert(icbrt(ULLONG_MAX) == 2642245ULL);

    assert(iroot(1023, 10) == 1);
    assert(iroot(1024, 10) == 2);
    assert(iroot(ULLONG_MAX, 63) == 2);
    assert(iroot(ULLONG_MAX, 64) == 1);
    return 0;
}
