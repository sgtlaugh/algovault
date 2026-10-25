/***
 *
 * Sparse Table
 * O(1) range queries on a static array, for operations where counting an element twice is harmless: min, max, gcd, and, or
 *
 * Complexity:
 *   - SparseTable: O(n log n) build time and memory, O(1) per query
 *   - LinearSparseTable: O(n) build time and memory, O(1) per query
 *
 * SparseTable<T> table(a): range minimum of a
 * SparseTable table(a, op): range op of a, for an op that ignores repeats like max or gcd
 * LinearSparseTable<T> table(a): range minimum of a
 * LinearSparseTable<T, greater<T>> table(a): range maximum of a
 * table.query(l, r): the answer over a[l..r], 0 <= l <= r < n
 *
 * Sum, xor and product count repeats, use disjoint_sparse_table.cpp for them
 * At n = 1e6 LinearSparseTable builds 3x faster in 12 MB instead of 72 MB, and its queries are about 2x slower
 *
***/

#include <bits/stdtr1c++.h>

using namespace std;

template <typename T>
struct MinOp{
    T operator()(const T& x, const T& y) const{
        return min(x, y);
    }
};

template <typename T, typename Op = MinOp<T>>
struct SparseTable{
    vector<vector<T>> dp;  /// dp[h][i] = op over a[i .. i + 2^h - 1]
    Op combine;

    SparseTable(const vector<T>& a, Op combine = Op()) : dp(1, a), combine(combine){
        int n = a.size();
        for (int h = 1; (1 << h) <= n; h++){
            dp.push_back(vector<T>(n - (1 << h) + 1));
            for (int i = 0; i + (1 << h) <= n; i++) dp[h][i] = combine(dp[h - 1][i], dp[h - 1][i + (1 << (h - 1))]);
        }
    }

    /// Two blocks of length 2^h cover [l, r], overlapping in the middle
    T query(int l, int r) const{
        assert(0 <= l && l <= r && r < (int)dp[0].size());
        int h = __lg(r - l + 1);
        return combine(dp[h][l], dp[h][r - (1 << h) + 1]);
    }
};

template <typename T, typename Compare = less<T>>
struct LinearSparseTable{
    static constexpr int B = 64;
    Compare compare;
    vector<T> val;
    vector<unsigned long long> mask;  /// bit k of mask[i]: val[i - k] beats every value after it up to i
    vector<vector<int>> table;        /// sparse table of block winners, by index

    LinearSparseTable(const vector<T>& a) : val(a), mask(a.size()){
        int n = val.size();
        unsigned long long stack = 0;  /// the monotonic stack of the last 64 positions, as a bitmask
        for (int i = 0; i < n; i++){
            stack <<= 1;
            /// the lowest set bit is the newest entry, pop while val[i] is at least as good
            while (stack && !compare(val[i - __builtin_ctzll(stack)], val[i])) stack &= stack - 1;
            stack |= 1;
            mask[i] = stack;
        }

        int blocks = (n + B - 1) / B;
        table.assign(1, vector<int>(blocks));
        for (int b = 0; b < blocks; b++) table[0][b] = best_ending_at(min(n - 1, b * B + B - 1), min(B, n - b * B));
        for (int k = 1; (1 << k) <= blocks; k++){
            table.push_back(vector<int>(blocks - (1 << k) + 1));
            for (int b = 0; b + (1 << k) <= blocks; b++) table[k][b] = better(table[k - 1][b], table[k - 1][b + (1 << (k - 1))]);
        }
    }

    T query(int l, int r) const{
        assert(0 <= l && l <= r && r < (int)val.size());
        if (r - l + 1 <= B) return val[best_ending_at(r, r - l + 1)];

        int res = better(best_ending_at(l + B - 1), best_ending_at(r));
        int x = l / B + 1, y = r / B - 1;
        if (x <= y){
            int k = __lg(y - x + 1);
            res = better(res, better(table[k][x], table[k][y - (1 << k) + 1]));
        }
        return val[res];
    }

    int better(int i, int j) const{
        return compare(val[j], val[i]) ? j : i;
    }

    /// Index of the best value in [r - len + 1, r], len <= 64: the oldest stack entry inside the window
    int best_ending_at(int r, int len = B) const{
        unsigned long long m = len == B ? mask[r] : mask[r] & ((1ULL << len) - 1);
        return r - __lg(m);
    }
};

template <typename T>
struct MaxOp{
    T operator()(const T& x, const T& y) const{
        return max(x, y);
    }
};

int main(){
    int n, m, t, d;

    scanf("%d", &n);
    vector<int> a(n);
    for (auto& v : a) scanf("%d", &v);

    vector<int> gap(max(1, n - 1), 0);
    for (int j = 0; j + 1 < n; j++) gap[j] = a[j + 1] - a[j];
    SparseTable<int, MaxOp<int>> table(gap);

    scanf("%d", &m);
    while (m--){
        scanf("%d %d", &t, &d);
        int k = upper_bound(a.begin(), a.end(), t) - a.begin() - 1;
        int lo = 0, hi = k;
        while (lo < hi){
            int mid = (lo + hi) / 2;
            if (table.query(mid, k - 1) <= d) hi = mid;
            else lo = mid + 1;
        }
        printf("%d\n", lo + 1);
    }
    return 0;
}
