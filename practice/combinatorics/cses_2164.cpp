/***
 *
 * Josephus Problem
 * n people numbered 1 to n stand in a circle, counting starts from 1 and every k-th person dies
 *
 * Complexity: O(min(m, k log n)), every step grows the circle by at least one person
 *
 * josephus(n, k, m): id in 1..n of the m-th person killed, 1 <= m <= n <= 4e18, 1 <= k <= 4e18
 *   josephus(n, k, n) is the survivor
 * For example if n = 10 and k = 3, then the people killed are 3, 6, 9, 2, 7, 1, 8, 5, 10, 4 respectively
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Runs the circle of size i from n - m + 1 up to n, jumping r sizes at once while no wrap-around happens
long long josephus(long long n, long long k, long long m){
    m = n - m;
    if (k <= 1) return n - m;

    for (long long i = m; i < n; ){
        long long r = i - m > k - 1 ? min((i - m + k - 2) / (k - 1), n - i) : 1;
        i += r;
        unsigned long long x = m + r * k;
        m = (x | i) >> 32 ? x % i : (unsigned)x % (unsigned)i;  /// 32-bit division when both fit, 1.5x faster for k >= n
    }
    return m + 1;
}

int main(){
    int q;
    scanf("%d", &q);
    while (q--){
        long long n, k;
        scanf("%lld %lld", &n, &k);
        printf("%lld\n", josephus(n, 2, k));
    }
    return 0;
}
