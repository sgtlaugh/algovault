/***
 *
 * Disjoint Sparse Table
 * O(1) range queries on a static array for any associative operation
 *
 * Complexity:
 *   - build: O(n log n) time and memory
 *   - query: O(1)
 *
 * DisjointSparseTable<T> table(a): range sums of a
 * DisjointSparseTable table(a, op): range folds with op, which must be associative but need not commute
 * table.query(l, r): op over a[l..r], 0 <= l <= r < n
 *
 * Use sparse_table.cpp for min, max or gcd, and a prefix array for plain sums
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Op = plus<T>>
struct DisjointSparseTable{
    vector<T> a;
    vector<vector<T>> dp;  /// dp[h][i]: op over a[i .. c - 1] if bit h of i is 0, over a[c .. i] if it is 1, c = middle of i's block of 2^(h + 1)
    Op combine;

    DisjointSparseTable(const vector<T>& a, Op combine = Op()) : a(a), combine(combine){
        int n = a.size();
        for (int h = 0; (1 << h) < n; h++){
            dp.push_back(a);
            for (int c = 1 << h; c < n; c += 2 << h){
                for (int i = c - 2; i >= c - (1 << h); i--) dp[h][i] = combine(a[i], dp[h][i + 1]);
                for (int i = c + 1; i < min(n, c + (1 << h)); i++) dp[h][i] = combine(dp[h][i - 1], a[i]);
            }
        }
    }

    /// l and r first differ at bit h, so l is in the left half of a block of 2^(h + 1) and r in its right half
    T query(int l, int r) const{
        assert(0 <= l && l <= r && r < (int)a.size());
        if (l == r) return a[l];
        int h = __lg(l ^ r);
        return combine(dp[h][l], dp[h][r]);
    }
};

struct Min{
    int operator()(int x, int y) const{
        return min(x, y);
    }
};

int main(){
    int n, q;
    scanf("%d %d", &n, &q);

    vector<int> a(n);
    for (auto& x: a) scanf("%d", &x);

    DisjointSparseTable<int, Min> table(a);
    while (q--){
        int l, r;
        scanf("%d %d", &l, &r);
        printf("%d\n", table.query(l - 1, r - 1));
    }

    return 0;
}
