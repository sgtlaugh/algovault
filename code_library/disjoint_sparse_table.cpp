/***
 *
 * Disjoint Sparse Table For RMQ on a static array
 * Idempotency is necessary for the regular sparse table
 * Because we answer queries in O(1) allowing partial overlaps which doesn't change the result
 * Some examples are finding minimum or maximum
 *
 * Disjoint Sparse Table can answer RMQ on a static array in O(1), even when they are not idempotent
 * Examples are sum, xor, product, etc (yes we can just use prefix sum for sum or prefix xor for xor)
 * Further reading - https://codeforces.com/blog/entry/87940
 *
 * Time and space complexity is the same as in regular sparse table
 * O(n log n) to build
 * O(1) to query
 *
 * Note, we can still do this in O(log n) using regular sparse tables
 * And then there are sqrt trees...
 * Disjoint sparse tables should be rarely useful in CP, nonetheless they are quite cool
 *
***/

#include <bits/stdc++.h>

using namespace std;

template <typename T>
struct DisjointST{
    vector <T> dp[32];
    vector <int> log_2;

    /// defined for xor by default, change as required
    /// only changing this function and using an appropriate identity element should suffice for all cases

    T combine(const T& x, const T& y){
        return x ^ y;
    }

    /// basically combine(x, identity_elem) = x for all x
    /// it is 0 for xor, 1 for product, -infinity for maximum and so on

    DisjointST(const vector<T> &ar, const T& identity_elem){
        int i, c, h, l, n, m;
        n = (int)ar.size(), m = 2 * n + 1;

        log_2.resize(m, 0);
        for (i = 2; i < m; i++) log_2[i] = log_2[i >> 1] + 1;

        for (h = 0; h <= log_2[n]; h++){
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

    T query(int l, int r){
        assert(l <= r);
        int h = log_2[l ^ (r + 1)];
        return combine(dp[h][l], dp[h][r + 1]);
    }
};

int main(){
    vector<int> v = {5, 6, 1, 13, 7, 4, 9, 66, 23};
    auto rmq = DisjointST<int>(v, 0);

    assert(rmq.query(0, 0) == 5);
    assert(rmq.query(0, 2) == 2);
    assert(rmq.query(1, 4) == 13);
    assert(rmq.query(3, 5) == 14);
    assert(rmq.query(4, 7) == 72);
    assert(rmq.query(2, 6) == 6);
    assert(rmq.query(1, 8) == 85);

    v.clear();
    mt19937 rng(0);

    const int n = 2000000;
    for (int i = 0; i < n; i++){
        v.push_back(rng() % 1000000000);
    }

    clock_t start = clock();
    rmq = DisjointST<int>(v, 0);
    fprintf(stderr, "\nTime taken to build = %0.6f\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Took 0.105 s locally
    return 0;
}
