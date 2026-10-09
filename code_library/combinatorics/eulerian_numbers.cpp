/***
 *
 * Eulerian Numbers - https://oeis.org/A008292
 * A(n, k) modulo m, the number of permutations of 1 to n in which exactly k elements are greater than their previous element
 *
 * Complexity: O(n^2) time and memory to build the triangle, O(1) per query
 *
 * EulerianNumbers eul(n, m) builds rows 0..n with the recurrence
 *     A(n, k) = (k + 1) A(n - 1, k) + (n - k) A(n - 1, k - 1), any modulus 1 <= m <= 2^31 - 1
 * eul.get(n, k) returns A(n, k) mod m, 0 when k < 0 or k >= n, and A(0, 0) = 1 by convention
 *
 * Eulerian triangle for n = 1 to 7 and k = 0 to n - 1 below
 *
 * 1
 * 1 1
 * 1 4 1
 * 1 11 11 1
 * 1 26 66 26 1
 * 1 57 302 302 57 1
 * 1 120 1191 2416 1191 120 1
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
    EulerianNumbers eul(1000);
    assert(eul.get(0, 0) == 1);
    assert(eul.get(1, 0) == 1 && eul.get(1, 1) == 0);
    assert(eul.get(4, 0) == 1 && eul.get(4, 1) == 11 && eul.get(4, 2) == 11 && eul.get(4, 3) == 1 && eul.get(4, 4) == 0);
    assert(eul.get(7, 2) == 1191 && eul.get(7, 3) == 2416 && eul.get(7, 6) == 1);
    assert(eul.get(5, -1) == 0 && eul.get(5, 9) == 0);
    assert(eul.get(1000, 500) == 948656644);

    EulerianNumbers big(1000, 2147483647);
    assert(big.get(1000, 1) == 2147482902 && big.get(1000, 500) == 958447772);

    EulerianNumbers small(7, 7);
    for (int k = 0; k < 7; k++) assert(small.get(7, k) == 1);
    assert(small.get(4, 1) == 4);

    EulerianNumbers unit(5, 1);
    assert(unit.get(0, 0) == 0 && unit.get(5, 2) == 0);

    return 0;
}
