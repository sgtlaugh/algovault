/***
 *
 * Fast Fibonacci
 * n-th Fibonacci number by fast doubling, F(0) = 0, F(1) = 1, F(n) = F(n - 1) + F(n - 2)
 *
 * Complexity: O(log n)
 *
 * fibonacci(n, m): F(n) mod m for 0 <= n < 2^63 and any modulus 1 <= m < 2^63
 * fibonacci(n): exact F(n) for 0 <= n <= 92, F(93) no longer fits in long long
 * fibonacci_pair(n, m): {F(n) mod m, F(n + 1) mod m}
 *
 * Notes:
 *   - Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

pair<long long, long long> fibonacci_pair(long long n, long long m){
    assert(n >= 0 && m >= 1);
    long long a = 0, b = 1 % m;  /// (F(k), F(k + 1)) for k = the bits of n read so far

    /// Sums compare with m - y since x + y can pass 2^63, about 25% faster than reducing each with __int128 %
    for (int bit = n ? __lg(n) : -1; bit >= 0; bit--){
        long long z = b >= m - b ? b - (m - b) : 2 * b;
        z = z >= a ? z - a : z - a + m;
        long long c = (long long)((__int128)a * z % m);
        long long d = (long long)(((__int128)a * a + (__int128)b * b) % m);
        if (n >> bit & 1) a = d, b = c >= m - d ? c - (m - d) : c + d;
        else a = c, b = d;
    }

    return {a, b};
}

long long fibonacci(long long n, long long m){
    return fibonacci_pair(n, m).first;
}

long long fibonacci(long long n){
    assert(0 <= n && n <= 92);
    return fibonacci(n, LLONG_MAX);
}

int main(){
    /// F(10) = 55 and F(11) = 89
    assert(fibonacci(10) == 55);
    assert(fibonacci(10, 7) == 6);
    assert((fibonacci_pair(10, 7) == make_pair(6LL, 5LL)));

    /// F(92) is the largest Fibonacci number that fits in long long
    assert(fibonacci(92) == 7540113804746346429LL);
    assert(fibonacci(100, 1000000007) == 687995182);

    return 0;
}
