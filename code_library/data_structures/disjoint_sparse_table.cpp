/***
 *
 * Disjoint Sparse Table
 * O(1) range queries on a static array for any associative operation: sum, xor, product mod p, matrix product
 *
 * Complexity: O(n log n) build and memory, O(1) per query
 *
 * DisjointSparseTable<T> table(a): table.query(l, r) is a[l] + ... + a[r], 0-based, l <= r
 * DisjointSparseTable table(a, op): the same with any associative op, it need not commute or have an identity
 *
 * Example, a = {5, 6, 1, 13, 7}:
 *     DisjointSparseTable<int> sum(a);
 *     sum.query(1, 3);     // 20
 *     sum.query(4, 4);     // 7
 *     DisjointSparseTable prod(a, [](int x, int y){ return x * y % 1000; });
 *     prod.query(0, 3);    // 390
 *
 * For min, max or gcd sparse_table.cpp is simpler and faster, for a plain sum a prefix array is enough
 * Further reading - https://codeforces.com/blog/entry/87940
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

int main(){
    /// The header example
    vector<int> a = {5, 6, 1, 13, 7};
    DisjointSparseTable<int> sum(a);
    assert(sum.query(1, 3) == 20);
    assert(sum.query(4, 4) == 7);
    DisjointSparseTable prod(a, [](int x, int y){ return x * y % 1000; });
    assert(prod.query(0, 3) == 390);

    /// Strings concatenate in order, so a wrong operand order would show
    DisjointSparseTable<string> words({"a", "b", "c", "d", "e", "f", "g", "h", "i"});
    assert(words.query(1, 3) == "bcd");
    assert(words.query(3, 8) == "defghi");
    assert(words.query(0, 8) == "abcdefghi");

    return 0;
}
