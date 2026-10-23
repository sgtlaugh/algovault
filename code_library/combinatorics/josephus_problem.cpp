/***
 *
 * Josephus Problem
 * n people numbered 1 to n stand in a circle, counting starts from 1 and every k-th person dies
 *
 * Complexity: O(min(m, k log n))
 *
 * josephus(n, k, m): id in 1..n of the m-th person killed, 1 <= m <= n <= 4e18, 1 <= k <= 4e18
 *
 * Notes:
 *   - josephus(n, k, n) is the survivor
 *   - n = 10 and k = 3 kill 3, 6, 9, 2, 7, 1, 8, 5, 10, 4 in that order
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

/// Runs the circle of size i from n - m + 1 up to n, jumping r sizes at once while no wrap-around happens
long long josephus(long long n, long long k, long long m){
    assert(1 <= m && m <= n && k >= 1);
    m = n - m;
    if (k <= 1) return n - m;

    for (long long i = m; i < n; ){
        long long r = i - m > k - 1 ? min((i - m + k - 2) / (k - 1), n - i) : 1;
        i += r;
        unsigned long long x = m + r * k;
        /// 32-bit division when both fit, 1.5x faster for k >= n
        m = (x | i) >> 32 ? x % i : (unsigned)x % (unsigned)i;
    }
    return m + 1;
}

int main(){
    vector<long long> order = {3, 6, 9, 2, 7, 1, 8, 5, 10, 4};
    for (int m = 1; m <= 10; m++) assert(josephus(10, 3, m) == order[m - 1]);

    assert(josephus(1, 5, 1) == 1);
    assert(josephus(7, 1, 4) == 4);
    assert(josephus(41, 3, 41) == 31);
    assert(josephus(100, 2, 100) == 73);  /// 100 = 64 + 36, the survivor for k = 2 is 2 * 36 + 1
    return 0;
}
