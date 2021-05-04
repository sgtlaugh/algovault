/***
 *
 * Disjoint Sparse Table For RMQ on a static array
 * Idempotency is necessary for the regular sparse table
 * Because we answer queries in O(1) allowing partial overlaps which doesn't change the result
 * Some examples are finding minimum or maximum
 *
 * Disjoint Sparse Table can answer RMQ on a static array in O(1), even when they are not idempotent
 * Examples are sum, xor, product (they should be associative)
 * Further reading - https://codeforces.com/blog/entry/87940
 *
 * Time and space complexity is the same as in regular sparse table
 * O(n log n) to build
 * O(1) to query
 *
***/

#include <bits/stdc++.h>

using namespace std;

/// defined for xor by default, change as required
struct Node{
    int value;

    /// neutral or identity element
    /// basically, x + i = x for all x (associative operation)
    /// where i is the identity element and + is the accumulator operator
    /// it is 0 for xor, 1 for product, -infinity for maximum and so on
    Node(): value(0) {}

    Node(long long value) : value(value) {}

    Node operator+ (const Node &other){ /// accumulator or merge operator
        return Node(value ^ other.value);
    }
};

struct DisjointST{
    vector <int> log_2;
    vector <Node> dp[32];

    template <typename T>
    DisjointST(const vector<T> &v){
        int i, c, h, l, n, m;
        n = (int)v.size(), m = 2 * n + 1;

        vector <Node> ar(n + 1);
        for (i = 0; i < n; i++) ar[i] = Node(v[i]);

        log_2.resize(m, 0);
        for (i = 2; i < m; i++) log_2[i] = log_2[i >> 1] + 1;

        for (h = 0; h <= log_2[n]; h++){
            dp[h].resize(n + 1);
            for (c = l = 1 << h; c < n + l; c += (l << 1)){
                for (i = c + 1; i < c + l && i <= n; i++){
                    dp[h][i] = dp[h][i - 1] + ar[i - 1];
                }

                for (i = min(n, c) - 1; i >= c - l; i--){
                    dp[h][i] = ar[i] + dp[h][i + 1];
                }
            }
        }
    }

    Node query(int l, int r){
        assert(l <= r);
        int h = log_2[l ^ (r + 1)];
        return dp[h][l] + dp[h][r + 1];
    }
};

int main(){
    vector<int> v = {5, 6, 1, 13, 7, 4, 9, 66, 23};
    auto rmq = DisjointST(v);

    assert(rmq.query(0, 0).value == 5);
    assert(rmq.query(0, 2).value == 2);
    assert(rmq.query(1, 4).value == 13);
    assert(rmq.query(3, 5).value == 14);
    assert(rmq.query(4, 7).value == 72);
    assert(rmq.query(2, 6).value == 6);
    assert(rmq.query(1, 8).value == 85);

    v.clear();
    mt19937 rng(0);

    const int n = 2000000;
    for (int i = 0; i < n; i++){
        v.push_back(rng() % 1000000000);
    }

    clock_t start = clock();
    rmq = DisjointST(v);
    fprintf(stderr, "\nTime taken to build = %0.6f\n", (clock() - start) / (double)CLOCKS_PER_SEC);  /// Took 0.105 s locally
    return 0;
}
