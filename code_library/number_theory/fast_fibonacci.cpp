/***
 *
 * Fast Fibonacci
 * n-th Fibonacci number by fast doubling, F(0) = 0, F(1) = 1, F(n) = F(n - 1) + F(n - 2)
 *
 * Complexity: O(log n)
 *
 * fibonacci(n, m): F(n) mod m for 0 <= n < 2^63 and any modulus 1 <= m < 2^63
 * fibonacci_pair(n, m): {F(n) mod m, F(n + 1) mod m}
 * fibonacci(n): exact F(n) for 0 <= n <= 92, F(93) no longer fits in long long
 *
 * F(2k) = F(k) (2 F(k + 1) - F(k)), F(2k + 1) = F(k)^2 + F(k + 1)^2, walking the bits of n from the top
 * Requires __int128 (64-bit GCC or Clang)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

pair<long long, long long> fibonacci_pair(long long n, long long m){
    assert(n >= 0 && m >= 1);
    auto mul = [&](long long a, long long b){ return (long long)((__int128)a * b % m); };
    long long a = 0, b = 1 % m;  /// (F(k), F(k + 1)) for k = the bits of n read so far

    for (int bit = n ? 63 - __builtin_clzll(n) : -1; bit >= 0; bit--){
        long long c = mul(a, ((__int128)2 * b - a + m) % m);
        long long d = (long long)(((__int128)a * a + (__int128)b * b) % m);
        if (n >> bit & 1) a = d, b = (long long)(((__int128)c + d) % m);  /// c + d can exceed 2^63 when m is close to it
        else a = c, b = d;
    }

    return {a, b};
}

long long fibonacci(long long n, long long m){
    return fibonacci_pair(n, m).first;
}

long long fibonacci(long long n){
    assert(0 <= n && n <= 92);
    if (n == 0) return 0;
    long long a = 0, b = 1;
    for (long long i = 1; i < n; i++) tie(a, b) = make_pair(b, a + b);  /// stops at F(n), F(93) would overflow
    return b;
}

int main(){
    const long long MOD = 1000000007;
    vector<long long> first = {0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55};

    for (int n = 0; n <= 10; n++) assert(fibonacci(n, MOD) == first[n] && fibonacci(n) == first[n]);
    assert(fibonacci(92) == 7540113804746346429LL);
    assert(fibonacci(92, LLONG_MAX) == 7540113804746346429LL);
    assert(fibonacci(100, MOD) == 687995182);
    assert(fibonacci(10, 1) == 0 && fibonacci(0, 1) == 0);

    assert((fibonacci_pair(10, 7) == make_pair(55LL % 7, 89LL % 7)));

    return 0;
}
