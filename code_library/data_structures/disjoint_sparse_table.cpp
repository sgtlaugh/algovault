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
 * Notes:
 *   - Use a prefix array for plain sums in O(n) memory, sparse_table.cpp for min, max or gcd
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T, typename Op = plus<T>>
struct DisjointSparseTable{
    vector<T> a;
    /// dp[h][i]: op over a[i .. c - 1] if bit h of i is 0, else over a[c .. i], c = middle of i's block of 2^(h + 1)
    vector<vector<T>> dp;
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

int main(){
    /// Sum and product over the same array
    vector<int> a = {5, 6, 1, 13, 7};
    DisjointSparseTable<int> sum(a);
    assert(sum.query(1, 3) == 20);
    assert(sum.query(4, 4) == 7);

    DisjointSparseTable prod(a, [](int x, int y){ return x * y % 1000; });
    assert(prod.query(0, 3) == 390);

    return 0;
}
