/***
 *
 * Faulhaber's Formula
 * Sum of k-th powers 1^k + 2^k + ... + n^k modulo a prime through Stirling numbers of the second kind
 *
 * Complexity: O(max_k^2) time and memory to build, O(k) per query
 *
 * Faulhaber<MOD> f(max_k) builds the Stirling triangle S(i, j) for i <= max_k
 *     MOD must be a prime > max_k + 1, the default is 10^9 + 7
 * f.sum(n, k) for 0 <= k <= max_k and any 0 <= n <= LLONG_MAX
 *     sum of i^k = sum over j of S(k, j) (n + 1) n (n - 1) ... (n + 1 - j) / (j + 1)
 *
 * Example:
 *   Faulhaber f(1000);
 *   f.sum(5, 2);  // 1 + 4 + 9 + 16 + 25 = 55
 *   Faulhaber<998244353> g(1000);
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <int MOD = 1000000007>
struct Faulhaber{
    int max_k;
    vector<vector<int>> stirling;
    vector<long long> inv;

    Faulhaber(int max_k) : max_k(max_k), stirling(max_k + 1), inv(max_k + 2){
        assert(max_k >= 0 && max_k + 1 < MOD);
        inv[1] = 1;
        for (int i = 2; i <= max_k + 1; i++) inv[i] = (MOD - MOD / i) * inv[MOD % i] % MOD;

        stirling[0] = {1};
        for (int i = 1; i <= max_k; i++){
            stirling[i].assign(i + 1, 0);
            for (int j = 1; j <= i; j++){
                long long keep = j < i ? (long long)stirling[i - 1][j] * j : 0;
                stirling[i][j] = (keep + stirling[i - 1][j - 1]) % MOD;
            }
        }
    }

    long long sum(long long n, int k) const{
        assert(n >= 0 && 0 <= k && k <= max_k);
        n %= MOD;
        /// S(0, 0) = 1 would count 0^0, so k = 0 is answered directly
        if (k == 0) return n;

        long long res = 0, falling = 1;
        for (int j = 0; j <= k; j++){
            falling = falling * ((n + 1 - j + MOD) % MOD) % MOD;
            res = (res + stirling[k][j] * falling % MOD * inv[j + 1]) % MOD;
        }
        return res;
    }
};

int main(){
    Faulhaber f(1000);
    assert(f.sum(9, 0) == 9);
    assert(f.sum(100, 1) == 5050);
    assert(f.sum(5, 2) == 55);     /// 1 + 4 + 9 + 16 + 25
    assert(f.sum(10, 3) == 3025);  /// 55^2, the cubes sum to the square of the first powers
    assert(f.sum(1000000007, 3) == 0);  /// (MOD (MOD + 1) / 2)^2 is a multiple of MOD
    assert(f.sum(1000000000000000000LL, 1000) == 486176152);

    Faulhaber<7> g(5);
    assert(g.sum(5, 2) == 6);      /// 55 mod 7
    assert(g.sum(10, 3) == 1);     /// 3025 mod 7
    assert(g.sum(8, 5) == 1);      /// 1^5 + ... + 6^5 = 0 mod 7, 7^5 = 0, 8^5 = 1
    return 0;
}
