/***
 *
 * Disjoint Sparse Table
 * O(1) range queries on a static array for any associative operation, idempotent or not
 *
 * Complexity: O(n log n) build and memory, O(1) per query, the same as the regular sparse table
 *
 * The regular sparse table answers with two overlapping blocks, so it only works for idempotent operations like min
 * This one splits every query into two disjoint precomputed halves, so sum, xor, product mod p and matrix products work
 * For sum or xor a prefix array is simpler, the table pays off for operations without an inverse
 * Further reading - https://codeforces.com/blog/entry/87940
 *
 * DisjointST<T, Op = plus<T>> table(values, identity, op), query(l, r) is Op over [l, r], 0-indexed and inclusive
 * identity must satisfy op(identity, x) = op(x, identity) = x: 0 for sum or xor, 1 for product, "" for concatenation
 * Op needs a const operator(), since query() is const, and need not be commutative
 *
 * Example:
 *   DisjointST<long long> sum(v);
 *   DisjointST prod(v, 1LL, [](long long a, long long b){ return a * b % MOD; });
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Op = plus<T>>
struct DisjointST{
    int n;
    vector<T> dp[32];
    Op combine;

    DisjointST(const vector<T> &ar, const T& identity_elem = T(), Op combine = Op()) : n(ar.size()), combine(combine){
        int i, c, h, l;

        for (h = 0; n && h <= __lg(n); h++){  /// __lg(0) is undefined
            dp[h].resize(n + 1, identity_elem);

            for (c = l = 1 << h; c < n + l; c += (l << 1)){
                for (i = c + 1; i <= min(n, c + l); i++){
                    dp[h][i] = combine(dp[h][i - 1], ar[i - 1]);
                }

                for (i = min(n, c) - 1; i >= c - l; i--){
                    dp[h][i] = combine(ar[i], dp[h][i + 1]);
                }
            }
        }
    }

    T query(int l, int r) const{
        assert(0 <= l && l <= r && r < n);
        int h = __lg(l ^ (r + 1));
        return combine(dp[h][l], dp[h][r + 1]);
    }
};

int main(){
    vector<int> v = {5, 6, 1, 13, 7, 4, 9, 66, 23};
    auto sum = DisjointST<int>(v, 0);

    assert(sum.query(0, 0) == 5);
    assert(sum.query(0, 2) == 12);
    assert(sum.query(1, 4) == 27);
    assert(sum.query(3, 5) == 24);
    assert(sum.query(4, 7) == 86);
    assert(sum.query(2, 6) == 34);
    assert(sum.query(1, 8) == 129);

    const long long MOD = 1000000007;
    DisjointST prod(vector<long long>{5, 6, 1, 13, 100000, 100000}, 1LL, [&](long long x, long long y){ return x * y % MOD; });
    assert(prod.query(0, 3) == 390 && prod.query(1, 2) == 6 && prod.query(4, 5) == 999999937 && prod.query(3, 3) == 13);

    v.clear();
    mt19937 rng(0);

    const int n = 2000000;
    for (int i = 0; i < n; i++){
        v.push_back(rng() % 1000);  /// partial sums span at most 2^20 elements, so they stay within int
    }

    clock_t start = clock();
    sum = DisjointST<int>(v, 0);
    fprintf(stderr, "\nTime taken to build = %0.6f\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Took 0.105 s locally

    return 0;
}
