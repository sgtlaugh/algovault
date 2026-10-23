/***
 *
 * Eulerian Numbers
 * A(n, k) modulo m, the number of permutations of 1..n with exactly k ascents
 *
 * Complexity: O(n^2) time and memory to build, O(1) per query
 *
 * EulerianNumbers eul(n, m): rows 0..n modulo any 1 <= m <= 2^31 - 1, m = 10^9 + 7 by default
 * eul.get(n, k): A(n, k) mod m, 0 when k < 0 or k >= n, and A(0, 0) = 1 by convention
 *
 * Notes:
 *   - An ascent is an element greater than the one before it
 *   - The triangle is OEIS A008292 (https://oeis.org/A008292)
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

struct EulerianNumbers{
    int mod;
    vector<vector<int>> dp;

    EulerianNumbers(int n, int mod = 1000000007): mod(mod), dp(n + 1){
        assert(n >= 0 && mod >= 1);
        for (int i = 0; i <= n; i++){
            dp[i].assign(i + 1, 0);
            dp[i][0] = 1 % mod;
            for (int k = 1; k < i; k++){
                dp[i][k] = ((long long)dp[i - 1][k] * (k + 1) + (long long)dp[i - 1][k - 1] * (i - k)) % mod;
            }
        }
    }

    int get(int n, int k) const{
        assert(0 <= n && n < (int)dp.size());
        if (k < 0 || k > n) return 0;
        return dp[n][k];
    }
};

int main(){
    /// Rows 0, 1, 4 and 7 of the triangle, and k outside [0, n) giving 0
    EulerianNumbers eul(1000);
    assert(eul.get(0, 0) == 1);
    assert(eul.get(1, 0) == 1);
    assert(eul.get(1, 1) == 0);
    assert(eul.get(4, 0) == 1);
    assert(eul.get(4, 1) == 11);
    assert(eul.get(4, 2) == 11);
    assert(eul.get(4, 3) == 1);
    assert(eul.get(4, 4) == 0);
    assert(eul.get(7, 2) == 1191);
    assert(eul.get(7, 3) == 2416);
    assert(eul.get(7, 6) == 1);
    assert(eul.get(5, -1) == 0);
    assert(eul.get(5, 9) == 0);
    assert(eul.get(1000, 500) == 948656644);

    EulerianNumbers big(1000, 2147483647);
    assert(big.get(1000, 1) == 2147482902);
    assert(big.get(1000, 500) == 958447772);

    EulerianNumbers small(7, 7);
    for (int k = 0; k < 7; k++) assert(small.get(7, k) == 1);
    assert(small.get(4, 1) == 4);

    EulerianNumbers unit(5, 1);
    assert(unit.get(0, 0) == 0);
    assert(unit.get(5, 2) == 0);

    return 0;
}
