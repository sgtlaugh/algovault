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
        int h = __lg(r - l);
        return combine(dp[h][l], dp[h][r - (1 << h) + 1]);
    }
};

int main(){
    vector<int> v = {5, 6, 1, 13, 7, 4, 9, 66, 23};
    auto rmq = SparseTable<int>(v);

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
