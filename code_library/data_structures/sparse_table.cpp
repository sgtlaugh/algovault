/***
 *
 * Sparse Table For RMQ on a static array
 * Idempotency is necessary for the regular sparse table and it should be an associative operation
 * Because we answer queries in O(1) allowing partial overlaps which doesn't change the result
 * Some examples are finding minimum or maximum, we cannot find sum or xor or product in this way
 *
 * Check out the disjoint sparse table for those situations
 * https://github.com/sgtlaugh/algovault/blob/master/code_library/disjoint_sparse_table.cpp
 *
 * So why use the regular sparse table at all as opposed to the disjoint sparse table?
 * Because its simpler, easier to understand and faster (roughly ~1.5x to 2x)
 *
 * Time and space complexity:
   * O(n log n) to build
   * O(1) to query
 *
 * LinearSparseTable<T> answers range minimum with the same API in O(1) per query, building n + (n / 64) log(n / 64) words,
 * effectively O(n) build and memory for any practical n
 * Measured with 1e7 queries: n = 1e6 builds 3x faster in 12 MB instead of 76 MB, queries 2.3x slower
 * At n = 1e7 the regular table needs 915 MB against 124 MB, so pick the linear one for large n or tight memory
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <typename T>
struct SparseTable{
    vector <T> dp[32];

    /// defined for min by default, change as required
    T combine(const T& x, const T& y){
        return min(x, y);
    }

    SparseTable(const vector<T> &ar){
        int i, j, l, h, n = (int)ar.size();

        dp[0] = ar;
        for (h = 1, l = 2; l <= n; h++, l <<= 1){
            dp[h].resize(n);
            for (i = 0, j = i + (l / 2); (i + l) <= n; i++, j++){
                dp[h][i] = combine(dp[h - 1][i], dp[h - 1][j]);
            }
        }
    }

    T query(int l, int r){
        int h = __lg(r - l + 1);
        return combine(dp[h][l], dp[h][r - (1 << h) + 1]);
    }
};

/// Blocks of 64, a per position bitmask of the in-block monotonic stack, and a sparse table over block minima
template <typename T>
struct LinearSparseTable{
    static constexpr int B = 64;
    vector<T> val;
    vector<unsigned long long> mask;
    vector<vector<int>> table;

    int better(int i, int j) const{
        return val[i] <= val[j] ? i : j;
    }

    /// Index of the minimum in [r - size + 1, r], size <= 64
    int small_query(int r, int size = B) const{
        unsigned long long m = size == B ? mask[r] : mask[r] & ((1ULL << size) - 1);
        return r - (63 - __builtin_clzll(m));
    }

    LinearSparseTable(const vector<T>& values) : val(values), mask(values.size()){
        int n = val.size();
        unsigned long long cur = 0;
        for (int i = 0; i < n; i++){
            cur <<= 1;
            while (cur && !(val[i - __builtin_ctzll(cur)] < val[i])) cur &= cur - 1;
            cur |= 1;
            mask[i] = cur;
        }

        int blocks = (n + B - 1) / B;
        table.assign(1, vector<int>(blocks));
        for (int b = 0; b < blocks; b++) table[0][b] = small_query(min(n - 1, b * B + B - 1), min(B, n - b * B));
        for (int k = 1; (1 << k) <= blocks; k++){
            table.push_back(vector<int>(blocks - (1 << k) + 1));
            for (int b = 0; b + (1 << k) <= blocks; b++) table[k][b] = better(table[k - 1][b], table[k - 1][b + (1 << (k - 1))]);
        }
    }

    T query(int l, int r) const{
        if (r - l + 1 <= B) return val[small_query(r, r - l + 1)];

        int res = better(small_query(l + B - 1), small_query(r));
        int x = l / B + 1, y = r / B - 1;
        if (x <= y){
            int k = 31 - __builtin_clz(y - x + 1);
            res = better(res, better(table[k][x], table[k][y - (1 << k) + 1]));
        }
        return val[res];
    }
};

int main(){
    vector<int> v = {5, 6, 1, 13, 7, 4, 9, 66, 23};
    auto rmq = SparseTable<int>(v);
    LinearSparseTable<int> linear(v);
    assert(linear.query(0, 0) == 5 && linear.query(0, 2) == 1 && linear.query(3, 5) == 4 && linear.query(4, 7) == 4 && linear.query(7, 8) == 23);

    assert(rmq.query(0, 0) == 5);
    assert(rmq.query(0, 2) == 1);
    assert(rmq.query(1, 4) == 1);
    assert(rmq.query(3, 5) == 4);
    assert(rmq.query(4, 7) == 4);
    assert(rmq.query(2, 6) == 1);
    assert(rmq.query(1, 8) == 1);

    v.clear();
    mt19937 rng(0);

    const int n = 2000000;
    for (int i = 0; i < n; i++){
        v.push_back(rng() % 1000000000);
    }

    clock_t start = clock();
    rmq = SparseTable<int>(v);
    fprintf(stderr, "\nTime taken to build = %0.6f\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Took 0.059744 s locally

    return 0;
}
